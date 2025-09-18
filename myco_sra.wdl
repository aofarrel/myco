version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.16.8/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.16.5/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.24/tasks/pull_fastqs.wdl" as sranwrp_pull
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.24/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/vcf_to_diff_wdl/process-metadata/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.3.0/tbprofiler_tasks.wdl" as profiler
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.3.0/theiagen_tbprofiler.wdl" as tbprofilerFQ_WF # fka earlyQC
import "https://raw.githubusercontent.com/aofarrel/goleft-wdl/0.1.3/goleft_functions.wdl" as goleft

# Copyright (C) 2025 Ash O'Farrell
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


workflow myco {
	input {
		File biosample_accessions

		Boolean just_like_2024                 = false
		Int     clean_average_q_score          = 29
		Boolean covstatsQC_skip_entirely       = true
		File?   mask_bedfile
		Boolean decontam_use_CDC_varpipe_ref   = false
		
		# QC stuff 
		Int     QC_max_pct_low_coverage_sites  =    20
		Int     QC_max_pct_unmapped            =     2
		Int     QC_min_mean_coverage           =    10
		Int     QC_min_q30                     =    90
		Boolean QC_soft_pct_mapped             = false
		Int     QC_this_is_low_coverage        =    10
		Int     quick_tasks_disk_size          =    10 
		Boolean guardrail_mode                 = true
		
		# shrink large samples
		Int     subsample_cutoff        =  450  # set to -1 to turn off subsampling entirely
		Int     subsample_seed          = 1965  # if you're trying replicate our results, leave this untouched!
	}

	parameter_meta {
		biosample_accessions: "File of BioSample accessions to pull, one accession per line"

		clean_average_q_score: "Trim reads with an average quality score below this value. Independent of QC_min_q30."
		covstatsQC_skip_entirely: "Should we skip covstats entirely?"
		mask_bedfile: "Bed file of regions to mask when making diff files (default: R00000039_repregions.bed)"

		QC_max_pct_low_coverage_sites: "Samples who have more than this percent (as int, 50 = 50%) of positions with coverage below QC_this_is_low_coverage will be discarded"
		QC_min_mean_coverage: "If covstats thinks MEAN coverage is below this, throw out this sample - not to be confused with TBProfiler MEDIAN coverage"
		QC_max_pct_unmapped: "If covstats thinks more than this percent of your sample (after decontam and cleaning) fails to map to H37Rv, throw out this sample."
		QC_min_q30: "Decontaminated samples with less than this percent (as int, 50 = 50%) of reads above qual score of 30 will be discarded."
		QC_soft_pct_mapped: "If true, a sample failing a percent mapped check (guardrail mode's TBProfiler check and/or covstats' check as per QC_max_pct_unmapped) will throw a non-fatal warning."
		QC_this_is_low_coverage: "Positions with coverage below this value will be masked in diff files"
		quick_tasks_disk_size: "Disk size in GB to use for quick file-processing tasks; increasing this might slightly speed up file localization"
		
		subsample_cutoff: "If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)"
		subsample_seed: "Seed used for subsampling with seqtk"
	}
	Boolean TBProf_on_bams_not_fastqs = just_like_2024
	String pass = "PASS" # used later... much later
	Float QC_max_pct_low_coverage_sites_float = QC_max_pct_low_coverage_sites / 100.0
	Int guardrail_subsample_cutoff = if guardrail_mode then 30000 else -1 # overridden by subsample_cutoff

	call sranwrp_processing.extract_accessions_from_file as get_sample_IDs {
		input:
			accessions_file = biosample_accessions,
			filter_na = true
	}

	scatter(biosample_accession in get_sample_IDs.accessions) {
		call sranwrp_pull.pull_fq_from_biosample as pull {
			input:
				biosample_accession = biosample_accession,
				fail_on_invalid = false,
				subsample_cutoff = select_first([subsample_cutoff, guardrail_subsample_cutoff]),
				subsample_seed = subsample_seed,
				tar_outputs = false
		} # output: pull.fastqs
		if(length(pull.fastqs)>1) {
    		Array[File] paired_fastqs=select_all(pull.fastqs)
  		}
	}

	call sranwrp_processing.cat_strings as merge_reports {
		input:
			strings = pull.results,
			out = "pull_reports.txt",
			disk_size = quick_tasks_disk_size
	}

	Array[Array[File]] pulled_fastqs = select_all(paired_fastqs)
	scatter(pulled_fastq in pulled_fastqs) {
		call clckwrk_combonation.clean_and_decontam_and_check as fastp_decontam_check {
			input:
				CDC_decontamination_reference = decontam_use_CDC_varpipe_ref,
				oldschool_docker = just_like_2024,
				unsorted_sam = true,
				reads_files = pulled_fastq,
				fastp_clean_avg_qual = clean_average_q_score,
				QC_min_q30 = QC_min_q30,
				preliminary_min_q30 = if guardrail_mode then 20 else 1,
				timeout_map_reads = if guardrail_mode then 120 else 0,
				timeout_decontam = if guardrail_mode then 300 else 0
				# no subsample cutoff here because that happens during the pull task
		}

		if(defined(fastp_decontam_check.decontaminated_fastq_1)) {
			# This region only executes if decontaminated fastqs exist. We can use this to coerce File? into File by using
			# select_first() where the first element is the File? we know must exist, and the second element is bogus.
    		File real_decontaminated_fastq_1=select_first([fastp_decontam_check.decontaminated_fastq_1, biosample_accessions])
    		File real_decontaminated_fastq_2=select_first([fastp_decontam_check.decontaminated_fastq_2, biosample_accessions])

			if(!(TBProf_on_bams_not_fastqs)) {
				call tbprofilerFQ_WF.ThiagenTBProfiler as theiagenTBprofilerFQ {
					input:
						fastq1 = real_decontaminated_fastq_1,
						fastq2 = real_decontaminated_fastq_2,
						soft_pct_mapped = QC_soft_pct_mapped,
						soft_depth = if guardrail_mode then false else true,
						minimum_depth = if guardrail_mode then 3 else 0,
						minimum_pct_mapped = if guardrail_mode then 10 else 0, # unlike covstats, this is a MINIMUM of % MAPPED
						sample = fastp_decontam_check.sample
				}
			}

			String tbprofiler_fq_status_or_bogus = select_first([theiagenTBprofilerFQ.status_code, "bogus"]) # prevent "cannot compare String? to String" error
			if(tbprofiler_fq_status_or_bogus == pass || TBProf_on_bams_not_fastqs) {
				call clckwrk_var_call.variant_call_one_sample_ref_included as variant_calling {
					input:
						reads_files = [real_decontaminated_fastq_1, real_decontaminated_fastq_2],
						tarball_bams_and_bais = false,
						timeout = if guardrail_mode then 600 else 0
				}
			}
		}
			
	}

	Array[File] minos_vcfs = flatten([select_all(variant_calling.adjudicated_vcf)])
	Array[File] final_bams = flatten([select_all(variant_calling.bam)])
	Array[File] final_bais = flatten([select_all(variant_calling.bai)])
	
	Array[Array[File]] bams_and_bais = [final_bams, final_bais]
	Array[Array[File]] bam_per_bai = transpose(bams_and_bais)
	
	scatter(vcfs_and_bams in zip(bam_per_bai, minos_vcfs)) {
	# scatter(vcfs_and_bams in zip(bam_per_bai, minos_vcfs)) is now sort of a three-way scatter:
	# * bam file accessible via vcfs_and_bams.left[0]
	# * bai file accessible via vcfs_and_bams.left[1]
	# * vcf file accessible via vcfs_and_bams.right
	
	# This relies on your WDL executor being consistent with how it orders arrays. That SHOULD always be the case per
	# the spec, but if things break catastrophically, let me save you some debug time: As of 2.9.2, clockwork-wdl's
	# ref-included version of the variant caller has an option to output the bams and bais as a tarball. You can use
	# that to recreate the simplier scatter of version 4.4.1 or earlier of myco. You will need to modify some tasks to
	# untar things, of course.
		if(!covstatsQC_skip_entirely) {
	
			# covstats to check coverage and percent mapped to reference
			call goleft.covstats as covstats {
				input:
					inputBamOrCram = vcfs_and_bams.left[0],
					allInputIndexes = [vcfs_and_bams.left[1]]
			}
			
			if((covstats.percentUnmapped < QC_max_pct_unmapped) || QC_soft_pct_mapped) {
				if(covstats.coverage > QC_min_mean_coverage) {
					
					# make diff files
					call diff.make_mask_and_diff as make_mask_and_diff_after_covstats {
						input:
							bam = vcfs_and_bams.left[0],
							vcf = vcfs_and_bams.right,
							min_coverage_per_site = QC_this_is_low_coverage,
							tbmf = mask_bedfile,
							max_ratio_low_coverage_sites_per_sample = QC_max_pct_low_coverage_sites_float
					}
				}
			}
		}
		
		if(covstatsQC_skip_entirely) {
		
			# make diff files
			call diff.make_mask_and_diff as make_mask_and_diff_no_covstats {
				input:
					bam = vcfs_and_bams.left[0],
					vcf = vcfs_and_bams.right,
					min_coverage_per_site = QC_this_is_low_coverage,
					tbmf = mask_bedfile,
					max_ratio_low_coverage_sites_per_sample = QC_max_pct_low_coverage_sites_float
			}
		}
		
		# TBProfiler (will run even if fails covstats qc)
		if(TBProf_on_bams_not_fastqs) {
			call profiler.tb_profiler_bam as profile_bam {
					input:
						bam = vcfs_and_bams.left[0]
			}
		}
	}

	# even though diffs and reports are technically optional outputs, this does work, and will avoid nulls in the final output
	Array[File] real_diffs = flatten([select_all(make_mask_and_diff_after_covstats.diff), select_all(make_mask_and_diff_no_covstats.diff)])
	Array[File] real_reports = flatten([select_all(make_mask_and_diff_after_covstats.report), select_all(make_mask_and_diff_no_covstats.report)])
	Array[File] real_masks = flatten([select_all(make_mask_and_diff_after_covstats.mask_file), select_all(make_mask_and_diff_no_covstats.mask_file)])

	# pull TBProfiler information, if we ran TBProfiler on bams
	
	# coerce optional types into required types (doesn't crash even if profile_bam didn't run)
	Array[String] coerced_bam_strains=select_all(profile_bam.sample_and_strain)
	Array[String] coerced_bam_resistances=select_all(profile_bam.sample_and_resistance)
	Array[String] coerced_bam_depths=select_all(profile_bam.sample_and_median_depth)
	
	# workaround for "defined(profile_bam.strain) is always true even if profile_bam didn't run" part of SOTHWO
	if(!(length(coerced_bam_strains) == 0)) {
	
		# if there is more than one sample, run some tasks to concatenate the outputs
		if(length(pulled_fastqs) != 1) {
			Array[String] bam_strains_with_header = flatten([["sample\tsublineage"], coerced_bam_strains])
			Array[String] bam_resista_with_header = flatten([["sample\tresistance"], coerced_bam_resistances])
			Array[String] bam_meddept_with_header = flatten([["sample\tmedn_depth"], coerced_bam_depths])
	
			call sranwrp_processing.cat_strings as collate_bam_strains {
				input:
					strings = bam_strains_with_header,
					out = "strain_reports.tsv",
					disk_size = quick_tasks_disk_size
			}
			
			call sranwrp_processing.cat_strings as collate_bam_resistance {
				input:
					strings = bam_resista_with_header,
					out = "resistance_reports.tsv",
					disk_size = quick_tasks_disk_size
			}
	
			call sranwrp_processing.cat_strings as collate_bam_depth {
				input:
					strings = bam_meddept_with_header,
					out = "depth_reports.tsv",
					disk_size = quick_tasks_disk_size
			}
		}
		
		# if there is only one sample, there's no need to run tasks
		if(length(pulled_fastqs) == 1) {
			String single_sample_tbprof_bam_depth      = coerced_bam_depths[0]
			String single_sample_tbprof_bam_resistance = coerced_bam_resistances[0]
			String single_sample_tbprof_bam_strain     = coerced_bam_strains[0]
		}
	}
  	
  	# pull TBProfiler information, if we ran TBProfiler on fastqs
  	
  	# coerce optional types into required types (doesn't crash if these are null)
	Array[String] coerced_fq_strains=select_all(theiagenTBprofilerFQ.sample_and_strain)
	Array[String] coerced_fq_resistances=select_all(theiagenTBprofilerFQ.sample_and_resistance)
	Array[String] coerced_fq_depths=select_all(theiagenTBprofilerFQ.sample_and_depth)
	
	# workaround for "defined(qc_fastq.strains) is always true" part of SOTHWO
	if(!(length(coerced_fq_strains) == 0)) {
	
		# if there is more than one sample, run some tasks to concatenate the outputs
		if(length(pulled_fastqs) != 1) {
			Array[String] fq_strains_with_header = flatten([["sample\tsublineage"], coerced_fq_strains])
			Array[String] fq_resista_with_header = flatten([["sample\tresistance"], coerced_fq_resistances])
			Array[String] fq_meddept_with_header = flatten([["sample\tmedn_depth"], coerced_fq_depths])

			call sranwrp_processing.cat_strings as collate_fq_strains {
				input:
					strings = fq_strains_with_header,
					out = "strain_reports.tsv",
					disk_size = quick_tasks_disk_size
			}
			
			call sranwrp_processing.cat_strings as collate_fq_resistance {
				input:
					strings = fq_resista_with_header,
					out = "resistance_reports.tsv",
					disk_size = quick_tasks_disk_size
			}
			
			call sranwrp_processing.cat_strings as collate_fq_depth {
				input:
					strings = fq_meddept_with_header,
					out = "depth_reports.tsv",
					disk_size = quick_tasks_disk_size
			}
		}
	
		# if there is only one sample, there's no need to run tasks
		if(length(pulled_fastqs) == 1) {
			String single_sample_tbprof_fq_depth      = coerced_fq_depths[0]
			String single_sample_tbprof_fq_resistance = coerced_fq_resistances[0]
			String single_sample_tbprof_fq_strain     = coerced_fq_strains[0]
		}
	}

	Array[String] columns = ["BioSample","raw_pct_above_q20","raw_pct_above_q30","raw_total_reads","post_cleaning_pct_above_q20","post_cleaning_pct_above_q30","post_decontam_pct_above_q20","post_decontam_pct_above_q30","post_decontam_total_reads","reads_is_contam","reads_TB_reference","reads_NTM","docker","status"]
	
	call sranwrp_processing.several_arrays_to_tsv as fastp_decont_report {
		input:
			row_keys = fastp_decontam_check.sample,
			column_keys = columns,
			value1 = fastp_decontam_check.q20_in,
			value2 = fastp_decontam_check.q30_in,
			value3 = fastp_decontam_check.reads_in,
			value4 = fastp_decontam_check.q20_postclean,
			value5 = fastp_decontam_check.q30_postclean,
			# cleaned_total_reads purposely excluded; it's borked
			value6 = fastp_decontam_check.q20_postdecon,
			value7 = fastp_decontam_check.q30_postdecon,
			value8 = fastp_decontam_check.reads_postdecon_per_fastp,
			value9 = fastp_decontam_check.reads_contam,
			value10 = fastp_decontam_check.reads_TB,
			value11 = fastp_decontam_check.reads_NTM,
			value12 = fastp_decontam_check.docker_used,
			value13 = fastp_decontam_check.error_code
	}
		
	output {
		File       download_report         = merge_reports.outfile
		File       fastp_decont_report_tsv = fastp_decont_report.tsv
		
		# raw files
		Array[File]  tbd_bais  = final_bais
		Array[File]  tbd_bams  = final_bams
		Array[File]  tbd_diffs = real_diffs
		Array[File]  tbd_masks = real_masks   # bedgraph
		Array[File]  tbd_vcfs  = minos_vcfs
		
		# metadata
		Array[File?] tbd_decontam_reports          = fastp_decontam_check.counts_out_tsv
		Array[File?] tbd_covstats_reports          = covstats.covstatsOutfile
		Array[File?] tbd_diff_reports              = real_reports
		Array[File?] tbd_tbprof_bam_jsons          = profile_bam.tbprofiler_json
		Array[File?] tbd_tbprof_bam_summaries      = profile_bam.tbprofiler_txt
		Array[File?] tbd_tbprof_fq_jsons           = theiagenTBprofilerFQ.tbprofiler_json
		Array[File?] tbd_tbprof_fq_looker          = theiagenTBprofilerFQ.tbprofiler_looker_csv
		Array[File?] tbd_tbprof_fq_laboratorian    = theiagenTBprofilerFQ.tbprofiler_laboratorian_report_csv
		Array[File?] tbd_tbprof_fq_lims            = theiagenTBprofilerFQ.tbprofiler_lims_report_csv
		
		# these outputs only exist if there are multiple samples
		File?        tbprof_bam_all_depths      = collate_bam_depth.outfile
		File?        tbprof_bam_all_strains     = collate_bam_strains.outfile
		File?        tbprof_bam_all_resistances = collate_bam_resistance.outfile
		File?        tbprof_fq_all_depths       = collate_fq_depth.outfile
		File?        tbprof_fq_all_strains      = collate_fq_strains.outfile
		File?        tbprof_fq_all_resistances  = collate_fq_resistance.outfile
		
		# these outputs only exist if we ran on a single sample
		String?      tbprof_bam_this_depth      = single_sample_tbprof_bam_depth
		String?      tbprof_bam_this_strain     = single_sample_tbprof_bam_strain
		String?      tbprof_bam_this_resistance = single_sample_tbprof_bam_resistance
		String?      tbprof_fq_this_depth       = single_sample_tbprof_fq_depth
		String?      tbprof_fq_this_strain      = single_sample_tbprof_fq_strain
		String?      tbprof_fq_this_resistance  = single_sample_tbprof_fq_resistance
	}
}

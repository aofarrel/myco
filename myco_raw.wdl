version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.16.8/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.16.5/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.24/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/vcf_to_diff_wdl/0.0.3/vcf_to_diff.wdl" as diff
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
		Array[Array[File]] paired_fastq_sets
		String date_pipeline_ran
		String? date_pipeline_previously_ran
		
		Boolean just_like_2024                 = false
		String? output_sample_name
		Boolean guardrail_mode                 = true
		Boolean low_resource_mode              = false
		Int     subsample_cutoff               =  -1 # note inconsistency with myco_sra!!
		
		Int     clean_average_q_score          = 29
		Boolean covstatsQC_skip_entirely       = true   # changed in myco 6.4.0
		Boolean decontam_use_CDC_varpipe_ref   = false  # changed in myco 6.3.0 -- # TODO: null op
		File?   mask_bedfile
		Int     QC_max_pct_low_coverage_sites  =    20
		Int     QC_max_pct_unmapped            =    10  # changed in myco 6.4.0
		Int     QC_min_mean_coverage           =    10  # CDC minimum: 50x
		Int     QC_min_q30                     =    80  # CDC minimum: 85%
		Boolean QC_soft_pct_mapped             = false
		Int     QC_this_is_low_coverage        =    10
		Int     quick_tasks_disk_size          =    10 
	}

	parameter_meta {
		clean_average_q_score: "Trim reads with an average quality score below this value. Independent of QC_min_q30."
		covstatsQC_skip_entirely: "Should we skip covstats entirely?"
		decontam_use_CDC_varpipe_ref: "If true, use CDC varpipe decontamination reference. If false, use CRyPTIC decontamination reference."
		guardrail_mode: "Implements about a half-dozen safeguards against extremely low-quality samples running for abnormally long times."
		mask_bedfile: "Bed file of regions to mask when making diff files (default: R00000039_repregions.bed)"
		output_sample_name: "Override all sample names with this string instead."
		paired_fastq_sets: "Nested array of paired fastqs, each inner array representing one samples worth of paired fastqs"
		QC_max_pct_low_coverage_sites: "Samples who have more than this percent (as int, 50 = 50%) of positions with coverage below QC_this_is_low_coverage will be discarded"
		QC_min_mean_coverage: "If covstats thinks MEAN coverage is below this, throw out this sample - not to be confused with TBProfiler MEDIAN coverage"
		QC_max_pct_unmapped: "If covstats thinks more than this percent of your sample (after decontam and cleaning) fails to map to H37Rv, throw out this sample."
		QC_min_q30: "Decontaminated samples with less than this percent (as int, 50 = 50%) of reads above qual score of 30 will be discarded."
		QC_soft_pct_mapped: "If true, a sample failing a percent mapped check (guardrail mode's TBProfiler check and/or covstats' check as per QC_max_pct_unmapped) will throw a non-fatal warning."
		QC_this_is_low_coverage: "Positions with coverage below this value will be masked in diff files"
		quick_tasks_disk_size: "Disk size in GB to use for quick file-processing tasks; increasing this might slightly speed up file localization"
	}
											  
	String pass = "PASS" # used later... much later
	Boolean tbprofiler_on_bam              = just_like_2024
	Int guardrail_subsample_cutoff = if guardrail_mode then 30000 else -1 # overridden by subsample_cutoff
	
	# flip some QC stuff around
	Float QC_max_pct_low_coverage_sites_float = QC_max_pct_low_coverage_sites / 100.0

	scatter(paired_fastqs in paired_fastq_sets) {
		call clckwrk_combonation.clean_and_decontam_and_check as decontam_each_sample {
			input:
				CDC_decontamination_reference = decontam_use_CDC_varpipe_ref,
				oldschool_docker = just_like_2024,
				unsorted_sam = true,
				force_rename_out = output_sample_name,
				reads_files = paired_fastqs,
				fastp_clean_avg_qual = clean_average_q_score,
				QC_min_q30 = QC_min_q30,
				preliminary_min_q30 = if guardrail_mode then 20 else 1,
				subsample_cutoff = select_first([subsample_cutoff, guardrail_subsample_cutoff]),
				timeout_map_reads = if guardrail_mode then 300 else 0,
				timeout_decontam = if guardrail_mode then 600 else 0,
				addldisk = if low_resource_mode then 10 else 100,
				memory = if low_resource_mode then 8 else 32
		}

		if(defined(decontam_each_sample.decontaminated_fastq_1)) {
			# This region only executes if decontaminated fastqs exist. We can use this to coerce File? into File by using
			# select_first() where the first element is the File? we know must exist, and the second element is bogus.
			# Originally I wanted to set the second element to something that isn't valid in hopes that may help catch 
			# otherwise silent errors should the behavior of Cromwell changes in the future with regard to defined()...
			# But it seems all members of a select_first() array are evaulated and checked for "is this valid" even if
			# the first (0 index) member of the select_first() array is defined (eg, is selected).
			File real_decontaminated_fastq_1=select_first([decontam_each_sample.decontaminated_fastq_1, paired_fastqs[0]])
			File real_decontaminated_fastq_2=select_first([decontam_each_sample.decontaminated_fastq_2, paired_fastqs[0]])
    		
			call tbprofilerFQ_WF.ThiagenTBProfiler as tbprofilerFQ {
				input:
					fastq1 = real_decontaminated_fastq_1,
					fastq2 = real_decontaminated_fastq_2,
					soft_pct_mapped = QC_soft_pct_mapped,
					soft_depth = if guardrail_mode then false else true,
					minimum_depth = if guardrail_mode then 3 else 0,
					minimum_pct_mapped = if guardrail_mode then 10 else 0, # unlike covstats, this is a MINIMUM of % MAPPED
					sample = decontam_each_sample.sample
			}
			# if this sample passes...
			if(tbprofilerFQ.status_code == pass) {
				call clckwrk_var_call.variant_call_one_sample_ref_included as variant_calling {
					input:
						reads_files = [real_decontaminated_fastq_1, real_decontaminated_fastq_2],
						tarball_bams_and_bais = false,
						timeout = if guardrail_mode then 600 else 0,
						addldisk = if low_resource_mode then 10 else 100,
						cpu = if low_resource_mode then 8 else 16,
						memory = if low_resource_mode then 8 else 32
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
	# bam file accessible via vcfs_and_bams.left[0]
	# bai file accessible via vcfs_and_bams.left[1]
	# vcf file accessible via vcfs_and_bams.right
	
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
		if(tbprofiler_on_bam) {
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
		if(length(paired_fastq_sets) != 1) {
	
			call sranwrp_processing.cat_strings as collate_bam_strains {
				input:
					strings = coerced_bam_strains,
					out = "strain_reports.txt",
					disk_size = quick_tasks_disk_size
			}
			
			call sranwrp_processing.cat_strings as collate_bam_resistance {
				input:
					strings = coerced_bam_resistances,
					out = "resistance_reports.txt",
					disk_size = quick_tasks_disk_size
			}
	
			call sranwrp_processing.cat_strings as collate_bam_depth {
				input:
					strings = coerced_bam_depths,
					out = "depth_reports.txt",
					disk_size = quick_tasks_disk_size
			}
		}
		
		# if there is only one sample, there's no need to run tasks
		if(length(paired_fastq_sets) == 1) {
			String single_sample_tbprof_bam_depth      = coerced_bam_depths[0]
			String single_sample_tbprof_bam_resistance = coerced_bam_resistances[0]
			String single_sample_tbprof_bam_strain     = coerced_bam_strains[0]
		}
	}
  	
  	# pull TBProfiler information, if we ran TBProfiler on fastqs
  	
  	# coerce optional types into required types (doesn't crash if these are null)
	Array[String] coerced_fq_strains=select_all(tbprofilerFQ.sample_and_strain)
	Array[String] coerced_fq_resistances=select_all(tbprofilerFQ.sample_and_resistance)
	Array[String] coerced_fq_depths=select_all(tbprofilerFQ.sample_and_depth)
	
	# workaround for "defined(qc_fastq.strains) is always true" part of SOTHWO
	if(!(length(coerced_fq_strains) == 0)) {
	
		# if there is more than one sample, run some tasks to concatenate the outputs
		if(length(paired_fastq_sets) != 1) {

			call sranwrp_processing.cat_strings as collate_fq_strains {
				input:
					strings = coerced_fq_strains,
					out = "strain_reports.txt",
					disk_size = quick_tasks_disk_size
			}
			
			call sranwrp_processing.cat_strings as collate_fq_resistance {
				input:
					strings = coerced_fq_resistances,
					out = "resistance_reports.txt",
					disk_size = quick_tasks_disk_size
			}
			
			call sranwrp_processing.cat_strings as collate_fq_depth {
				input:
					strings = coerced_fq_depths,
					out = "depth_reports.txt",
					disk_size = quick_tasks_disk_size
			}
		}
	
		# if there is only one sample, there's no need to run tasks
		if(length(paired_fastq_sets) == 1) {
			String single_sample_tbprof_fq_depth      = coerced_fq_depths[0]
			String single_sample_tbprof_fq_resistance = coerced_fq_resistances[0]
			String single_sample_tbprof_fq_strain     = coerced_fq_strains[0]
		}
	}
	
	#########################################
	# error reporting for Terra data tables #
	#########################################
	#
	# When running on a Terra data table, one instance of the workflow is created for every sample. This is in contrast to how
	# running one instance of the workflow to handle multiple samples. In the one-instance case, we can return an error code
	# for an individual sample as workflow-level output, which gets written to the Terra data table.
	
	# is there only one sample?
	if(length(paired_fastq_sets) == 1) {

		# did the decontamination step actually run? (note that defined() is not a robust check, but since this is the first task
		# in the workflow this should be okay for now)
		if(defined(decontam_each_sample.error_code)) {
		
			# Sidenote: If you have a task taking in decontam_each_sample.error_code, even after this defined check, it will fail
			# command instantiation with "Cannot interpolate Array[String?] into a command string with attribute set 
			# [PlaceholderAttributeSet(None,None,None,Some( ))]".
		
			# get the first (0th) value, eg only value since there's just one sample, and coerce it into type String
			String coerced_decontam_errorcode = select_first([decontam_each_sample.error_code[0], "WORKFLOW_ERROR_1_REPORT_TO_DEV"])
			
			# did the decontamination step return an error?
			if(!(coerced_decontam_errorcode == pass)) {          
				String decontam_ERR = coerced_decontam_errorcode
			}
		}
		
		# handle earlyQC (if it ran at all)
		
		Array[String] earlyQC_array_coerced = select_all(tbprofilerFQ.status_code)
		Array[String] earlyQC_errorcode_array = flatten([earlyQC_array_coerced, ["PASS"]]) # will fall back to PASS if earlyQC was skipped
		if(!(earlyQC_errorcode_array[0] == pass)) {          
			String earlyQC_ERR = earlyQC_errorcode_array[0]
		}

		# The WDL 1.0 spec does not say what happens if you give select_all() an array that only has optional values, but
		# the WDL 1.1 spec says you get an empty array. Thankfully, Cromwell handles 1.0-select_all() like the 1.1 spec.
		Array[String] errorcode_if_earlyQC = select_all(variant_calling.errorcode)
		
		# if the variant caller did not run, the fallback pass will be selected, even though the sample shouldn't be considered a pass, so
		# the final-final-final error code needs to have decontam's error come before the variant caller error.
		Array[String] varcall_errorcode_array = flatten([errorcode_if_earlyQC, ["PASS"]])
		if(!(varcall_errorcode_array[0] == pass)) {          
				String varcall_ERR = varcall_errorcode_array[0]
		}
		
		# handle covstats
		if (!covstatsQC_skip_entirely) {
			if (varcall_errorcode_array[0] == "PASS") {
				if(length(covstats.percentUnmapped) > 0) {
					# cannot use defined(covstats.percentUnmapped) as it is always true, so we instead
					# check if there are more than zero values in the output array for covstats
					Array[Float] percentsUnmapped = select_all(covstats.percentUnmapped)
					Float        percentUnmapped = percentsUnmapped[0]
					Array[Float] meanCoverages = select_all(covstats.coverage)
					Float        meanCoverage = meanCoverages[0]
					
					if((percentUnmapped > QC_max_pct_unmapped) && !(QC_soft_pct_mapped)) { 
						String too_many_unmapped = "COVSTATS_${percentUnmapped}_UNMAPPED_(MAX_${QC_max_pct_unmapped})"
						if(meanCoverage < QC_min_mean_coverage) {
							String double_bad = "COVSTATS_BOTH_${percentUnmapped}_UNMAPPED_(MAX_${QC_max_pct_unmapped})_AND_${meanCoverage}_MEAN_COVERAGE_(MIN_${QC_min_mean_coverage})"
						} 
					}
					if(meanCoverage < QC_min_mean_coverage) {
						String too_low_coverage = "COVSTATS_${meanCoverage}_MEAN_COVERAGE_(MIN_${QC_min_mean_coverage})"
					}
				}
			}
			String coerced_covstats_error = select_first([double_bad, too_low_coverage, too_many_unmapped, "PASS"])
			if(!(coerced_covstats_error == pass)) {          
					String covstats_ERR = coerced_covstats_error
			}
		}
		
		# handle vcf to diff
		# will use the same workaround as the variant caller
		Array[String] vcfdiff_errorcode_if_covstats = select_all(make_mask_and_diff_after_covstats.errorcode)
		Array[String] vcfdiff_errorcode_if_no_covstats = select_all(make_mask_and_diff_no_covstats.errorcode)
		Array[String] vcfdiff_errorcode_array = flatten([vcfdiff_errorcode_if_covstats, vcfdiff_errorcode_if_no_covstats, ["PASS"]])
		if(!(vcfdiff_errorcode_array[0] == pass)) {          
				String vcfdiff_ERR = vcfdiff_errorcode_array[0]
		}
		
		# final-final-final error code
		# earlyQC is at the end (but before PASS) to account for earlyQC_skip_QC = true
		String finalcode = select_first([decontam_ERR, varcall_ERR, covstats_ERR, vcfdiff_ERR, earlyQC_ERR, pass])
	}
	
	# miniwdl check will allow using just one flatten() here, but womtool will not. per the spec, flatten() isn't recursive.
	# TODO: this is still breaking in Cromwell!
	# Failed to evaluate 'warnings' (reason 1 of 1): Evaluating flatten(flatten([[select_all(tbprofilerFQ.warning_codes)], 
	# [select_all(warning_decontam)]])) failed: No coercion defined from wom value(s) '
	# [["EARLYQC_88.112_PCT_ABOVE_Q30_(MIN_0.9)", "EARLYQC_99.61_PCT_MAPPED_(MIN_99.995)"]]' of type 'Array[Array[String]]' to 'Array[String]'.
	#Array[String] warnings = flatten(flatten([[select_all(tbprofilerFQ.warning_codes)], [select_all(warning_decontam)]]))
		
	output {
		String tbd_status = select_first([finalcode, pass])
		String tbd_pipeline_run = select_first([date_pipeline_previously_ran, date_pipeline_ran])

		# decon/fastp metadata pulled out directly
		Float tbd_qc_q20_in = decontam_each_sample.q20_in[0]
		Float tbd_qc_q30_in = decontam_each_sample.q30_in[0]
		Int tbd_qc_reads_in = decontam_each_sample.reads_in[0]
		Int tbd_qc_mean_r1_len_in = decontam_each_sample.mean_r1_len_in[0]
		Int tbd_qc_mean_r2_len_in = decontam_each_sample.mean_r2_len_in[0]
		Float tbd_qc_q20_postclean = decontam_each_sample.q20_postclean[0]
		Float tbd_qc_q30_postclean = decontam_each_sample.q30_postclean[0]
		Int tbd_qc_reads_postclean_per_fastp = decontam_each_sample.reads_postclean_per_fastp[0]
		Int tbd_qc_mean_r1_len_postclean = decontam_each_sample.mean_r1_len_postclean[0]
		Int tbd_qc_mean_r2_len_postclean = decontam_each_sample.mean_r2_len_postclean[0]
		Float tbd_qc_pct_loss_cleaning_per_fastp = decontam_each_sample.pct_loss_cleaning_per_fastp[0]
		Int tbd_qc_reads_postclean_per_decon = decontam_each_sample.reads_postclean_per_decon[0]
		Int tbd_qc_reads_postdecon_per_decon = decontam_each_sample.reads_postdecon_per_decon[0]
		Int tbd_qc_reads_TB = decontam_each_sample.reads_TB[0]
		Int tbd_qc_reads_NTM = decontam_each_sample.reads_NTM[0]
		Int tbd_qc_reads_human = decontam_each_sample.reads_human[0]
		Int tbd_qc_reads_contam = decontam_each_sample.reads_contam[0]
		Float tbd_qc_pct_reads_TB_predecon = decontam_each_sample.pct_reads_TB_predecon[0]
		Float tbd_qc_pct_reads_NTM = decontam_each_sample.pct_reads_NTM[0]
		Float tbd_qc_pct_reads_human = decontam_each_sample.pct_reads_human[0]
		Float tbd_qc_pct_reads_TB_postdecon = decontam_each_sample.pct_reads_TB_postdecon[0]
		Float tbd_qc_pct_loss_decon_per_decon = decontam_each_sample.pct_loss_decon_per_decon[0]
		Float tbd_qc_pct_loss_total = decontam_each_sample.pct_loss_total[0]
		Float tbd_qc_pct_loss_decon_per_fastp = decontam_each_sample.pct_loss_decon_per_fastp[0]
		Float tbd_qc_q20_postdecon = decontam_each_sample.q20_postdecon[0]
		Float tbd_qc_q30_postdecon = decontam_each_sample.q30_postdecon[0]
		Int tbd_qc_reads_postdecon_per_fastp = decontam_each_sample.reads_postdecon_per_fastp[0]
		Int tbd_qc_mean_r1_len_postdecon = decontam_each_sample.mean_r1_len_postdecon[0]
		Int tbd_qc_mean_r2_len_postdecon = decontam_each_sample.mean_r2_len_postdecon[0]
		Float tbd_qc_duplication_rate = decontam_each_sample.duplication_rate[0]
		Int tbd_qc_reads_adapter_trimmed = decontam_each_sample.reads_adapter_trimmed[0]
		
		# theiagen!TBProfiler metadata pulled out directly
		Float? tbd_qc_median_depth_per_tbprof = tbprofilerFQ.median_depth[0]
		Float? tbd_qc_avg_depth_per_tbprof = tbprofilerFQ.avg_depth[0]
		Float? tbd_qc_pct_mapped_per_tbprof = tbprofilerFQ.pct_reads_mapped[0]
		Float? tbd_qc_pct_genome_covered = tbprofilerFQ.pct_genome_covered[0]
		Int?    tbd_n_dr_variants = tbprofilerFQ.n_dr_variants[0]
		Int?    tbd_n_other_variants = tbprofilerFQ.n_other_variants[0]
		String? tbd_resistance = tbprofilerFQ.resistance[0]
		String? tbd_strain_per_tbprof = tbprofilerFQ.strain[0]
		
		# covstats
		Float? tbd_qc_pct_unmapped_covstats = percentUnmapped
		
		# intermediate files
		Array[File]  tbd_bais  = final_bais
		Array[File]  tbd_bams  = final_bams
		Array[File]  tbd_diffs = real_diffs
		Array[File]  tbd_masks = real_masks   # bedgraph
		Array[File]  tbd_vcfs  = minos_vcfs
		
		# metadata files
		Array[File?] tbd_decontam_reports          = decontam_each_sample.counts_out_tsv
		Array[File?] tbd_covstats_reports          = covstats.covstatsOutfile
		Array[File?] tbd_diff_reports              = real_reports
		Array[File?] tbd_tbprof_bam_jsons          = profile_bam.tbprofiler_json
		Array[File?] tbd_tbprof_bam_summaries      = profile_bam.tbprofiler_txt
		Array[File?] tbd_tbprof_fq_jsons           = tbprofilerFQ.tbprofiler_json
		Array[File?] tbd_tbprof_fq_looker          = tbprofilerFQ.tbprofiler_looker_csv
		Array[File?] tbd_tbprof_fq_laboratorian    = tbprofilerFQ.tbprofiler_laboratorian_report_csv
		Array[File?] tbd_tbprof_fq_lims            = tbprofilerFQ.tbprofiler_lims_report_csv

		# useful debugging/run information (only valid iff this ran on only one sample)
		#Array[String] pass_or_warnings = if (length(warnings) > 0) then warnings else ["PASS"]
		#String? tbd_debug_decontam_ERR  = decontam_ERR
		#String? tbd_debug_earlyQC_ERR   = earlyQC_ERR
		#String? tbd_debug_varcall_ERR   = varcall_ERR
		#String? tbd_debug_covstats_ERR  = covstats_ERR
		#String? tbd_debug_vcfdiff_ERR   = vcfdiff_ERR
		#Array[String]? tbd_debug_vcfdiff_errorcode_if_covstats    = vcfdiff_errorcode_if_covstats
		#Array[String]? tbd_debug_vcfdiff_errorcode_if_no_covstats = vcfdiff_errorcode_if_no_covstats
		#Array[String]? tbd_debug_vcfdiff_errorcode_array          = vcfdiff_errorcode_array
		String tbd_clockwork_docker       = decontam_each_sample.docker_used[0]
	}
}


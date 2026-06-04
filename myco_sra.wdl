version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.16.10/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.16.10/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.2.2/tasks/pull_fastqs.wdl" as sranwrp_pull
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.2.2/tasks/processing_tasks.wdl" as sranwrp_processing
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
		File? biosample_accessions_file
		String biosample_accession_str   # if using biosample_accessions_file this can be an empty string

		Boolean just_like_2024                 = false
		Int     clean_average_q_score          = 29
		Boolean skip_covstats                  = true
		Boolean generate_download_report_file  = true
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
		biosample_accessions_file: "File of multiple BioSample accessions to pull, one accession per line. Recommended for non-Terra users. Overrides biosample_accession_str."
		biosample_accession_str: "String of one (1) BioSample accession to pull. Recommended for Terra users. If biosample_accessions_file exists, biosample_accession_str will be ignored."

		clean_average_q_score: "Trim reads with an average quality score below this value. Independent of QC_min_q30."
		skip_covstats: "Should we skip covstats entirely?"
		generate_download_report_file: "Generate file reporting all pulls (recommended if multi-sample batch that uses biosample_accessions_file)"
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

	if (defined(biosample_accessions_file)) {
		call sranwrp_processing.extract_accessions_from_file_with_fake_optional_input as get_sample_IDs {
			input:
				accessions_file = biosample_accessions_file,
				filter_na = true
		}
	}
	
	Array[String] biosample_accessions = select_first([get_sample_IDs.accessions, [biosample_accession_str]])
	
	scatter(biosample_accession in biosample_accessions) {
		call sranwrp_pull.pull_fq_from_biosample as pull {
			input:
				biosample_accession = biosample_accession,
				fail_on_invalid = false,
				subsample_cutoff = select_first([subsample_cutoff, guardrail_subsample_cutoff]),
				subsample_seed = subsample_seed,
				tar_outputs = false
		}
		if(length(pull.fastqs)>1) {
			Array[File] paired_fastqs=select_all(pull.fastqs)
		}
		if(!(length(pull.fastqs)>1)) {
			String pull_ERR = "ERROR_PULLING_FROM_SRA_SEE_DOWNLOAD_REPORT"
		}
	}

	if (generate_download_report_file) {
		call sranwrp_processing.cat_strings as merge_reports {
			input:
				strings = pull.results,
				out = "pull_reports.txt",
				disk_size = quick_tasks_disk_size
		}
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
				strip_all_underscores = true,
				preliminary_min_q30 = if guardrail_mode then 20 else 1,
				timeout_map_reads = if guardrail_mode then 120 else 0,
				timeout_decontam = if guardrail_mode then 300 else 0
				# no subsample cutoff here because that happens during the pull task
		}

		if(defined(fastp_decontam_check.decontaminated_fastq_1)) {
			# This region only executes if decontaminated fastqs exist. We can use this to coerce File? into File by using
			# select_first() where the first element is the File? we know must exist, and the second element is bogus.
			#
			# NOTE: biosample_accessions used to be an always-defined file, now it is an always-defined Array[String], so this
			# will warn FileCoercion in miniwdl check. However, it passes womtool, and as stated above the fallback is never
			# actually used, so we currently can get away with this. But if Cromwell's type checker ever becomes even stricter,
			# this might one day throw an error. If that happens, the easiest workaround is to force the user to input some
			# kind of random bogus file and use that as the fallback here and elsewhere, or to take the approach in myco_raw,
			# which in some ways is even more cursed.
			#
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

	# Will be used to determine if we need to concatenate QC information into a file or just present as-is
	Boolean is_single_sample_run = if (length(pulled_fastqs) == 1) then true else false

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
		if(!skip_covstats) {
	
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
		
		if(skip_covstats) {
		
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

	#########################################
	#      TBProfiler metadata handling     #
	#########################################
	# This next section follows this line of logic:
	# 1. Coerce TBProfiler's optional outputs into required types
	# ---> This prevents certain WDL crashes and, somehow, doesn't crash even if the files don't exist
	#      at runtime.
	# 2. Determine if we ran TBProfiler on bams/FQs
	# ---> Ideally we would know if TBProfiler ran on bams (which as a task is called profile_bam) by checking
	#      if one of the outputs of profile_bam is defined(). However, for some bloody reason (see my comments
	#      elsewhere about "SOTHWO", even if profile_bam didn't run, defined(profile_bam.strain) is always true!
	#      This seems to be a result of Cromwell creating empty arrays for all possible outputs within a scatter(),
	#      so the array itself is defined, it just has no values. Now, WHY Cromwell would premake empty arrays for
	#      tasks that will never run is beyond me... but it does (at least it did in 2023; even if this behavior has
	#      since changed that would mean I couldn't rely upon defined() not chaging b/n versions.)
	# ---> As a workaround, we check if the metadata we coerced in #1 is an array of a non-zero length.
	# ---> Can't we just check if length(profile_bam.sample_and_strain), ie the non-coerced version, has a non-zero
	#      length? Maybe. IIRC there is a miniwdl/Cromwell inconsistency so don't do that unless something breaks.
	# 3. Determine if we are running on one sample or multiple samples
	# ---> If we are running on multiple samples it is worth our time concatenating a bunch of lists into one
	#      metadata file. If we are running on just one sample this is not worth the compute cost/time.
	
	# coerce optional types into required types (doesn't crash even if profile_bam didn't run)
	Array[String] coerced_bam_strains=select_all(profile_bam.sample_and_strain)
	Array[String] coerced_bam_resistances=select_all(profile_bam.sample_and_resistance)
	Array[String] coerced_bam_depths=select_all(profile_bam.sample_and_median_depth)
	
	# workaround for "defined(profile_bam.strain) is always true even if profile_bam didn't run" part of SOTHWO
	if(!(length(coerced_bam_strains) == 0)) {
	
		# if there is more than one sample, run some tasks to concatenate the outputs
		if(!(is_single_sample_run)) {
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
		if(is_single_sample_run) {
			# these are tab-delimited and break Terra data tables, so they're no longer workflow outputs
			# TODO: either remove tbprof on bam or output its stats to the terra data table too
			String single_sample_tbprof_bam_samp_tab_depth      = coerced_bam_depths[0]
			String single_sample_tbprof_bam_samp_tab_resistance = coerced_bam_resistances[0]
			String single_sample_tbprof_bam_samp_tab_strain     = coerced_bam_strains[0]
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
		if(is_single_sample_run) {
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
	
		# if there is only one sample, there's no need to run tasks, just output directly
		if(is_single_sample_run) {
			# these are tab-delimited and break Terra data tables, so they're no longer workflow outputs
			String single_sample_tbprof_fq_samp_tab_depth      = coerced_fq_depths[0]
			String single_sample_tbprof_fq_samp_tab_resistance = coerced_fq_resistances[0]
			String single_sample_tbprof_fq_samp_tab_strain     = coerced_fq_strains[0]

			# these ones are not tab-delimited (nor are they coerced, but I think that's okay)
			Float?  single_sample_tbprof_fq_median_depth      = theiagenTBprofilerFQ.median_depth[0]
			Float?  single_sample_tbprof_fq_avg_depth         = theiagenTBprofilerFQ.avg_depth[0]
			Float?  single_sample_tbprof_fq_pct_mapped        = theiagenTBprofilerFQ.pct_reads_mapped[0]
			Float?  single_sample_tbprof_fq_pct_covered       = theiagenTBprofilerFQ.pct_genome_covered[0]
			Int?    single_sample_tbprof_fq_n_dr_variants     = theiagenTBprofilerFQ.n_dr_variants[0]
			Int?    single_sample_tbprof_fq_n_other_variants  = theiagenTBprofilerFQ.n_other_variants[0]
			String? single_sample_tbprof_fq_resistance        = theiagenTBprofilerFQ.resistance[0]
			String? single_sample_tbprof_fq_strain            = theiagenTBprofilerFQ.strain[0]

		}
	}

	Array[String] columns = ["BioSample","raw_pct_above_q20","raw_pct_above_q30","raw_total_reads","post_cleaning_pct_above_q20","post_cleaning_pct_above_q30","post_decontam_pct_above_q20","post_decontam_pct_above_q30","post_decontam_total_reads","reads_is_contam","reads_TB_reference","reads_NTM","docker","status"]

	if (!(is_single_sample_run)) {
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
	}
	
	

	#########################################
	# error reporting for Terra data tables #
	#########################################
	#
	# When running on a Terra data table, one instance of the workflow is created for every sample. (As opposed to running one
	# one instance of the workflow that handles multiple samples, which is what we did for NCBI SRA samples.) In the one-instance
	# case, we can return an error code for an individual sample as workflow-level output, which gets written to the Terra data table.
	
	# is there only one sample?
	if(is_single_sample_run) {

		String zeroth_sample_pull_code = pull.results[0]

		# did the decontamination step actually run? (note that defined() is not a robust check, but since this is the first task
		# in the workflow this should be okay for now)
		if(defined(fastp_decontam_check.error_code)) {
		
			# Sidenote: If you have a task taking in fastp_decontam_check.error_code, even after this defined check, it will fail
			# command instantiation with "Cannot interpolate Array[String?] into a command string with attribute set 
			# [PlaceholderAttributeSet(None,None,None,Some( ))]".
		
			# get the first (0th) value, eg only value since there's just one sample, and coerce it into type String
			String coerced_decontam_errorcode = select_first([fastp_decontam_check.error_code[0], "WORKFLOW_ERROR_1_REPORT_TO_DEV"])
			
			# did the decontamination step return an error?
			if(!(coerced_decontam_errorcode == pass)) {          
				String decontam_ERR = coerced_decontam_errorcode
			}
		}
		
		# handle Theiagen TBPRofiler, which serves as earlyQC (if it ran at all)
		
		Array[String] earlyQC_array_coerced = select_all(theiagenTBprofilerFQ.status_code)
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
		if (!skip_covstats) {
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
		# because skipping FQ TBProfiler via earlyQC_skip_QC is no longer an option, its code is no longer at the end
		# because covstats_ERR is undefined if !skip_covstats, covstats_ERR should not short-circuit to pass
		# TODO: because zeroth_sample_pull_code is defined regardless of pass/fail, if it's at the front we will never
		# get error codes and if it's before pass we will never fall back to pass
		# --> Adding pull_ERR throws "unable to unify array types," unsure why, maybe it needs to be coerced non-optional?
		String finalcode = select_first([decontam_ERR, earlyQC_ERR, varcall_ERR, covstats_ERR, vcfdiff_ERR, pass])
	}
	String multi_sample_status_code = "multi-sample run"
		
	output {
		File?      download_report         = merge_reports.outfile
		String     tbd_status              = select_first([finalcode, multi_sample_status_code])
		
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
		File?        fastp_decont_report_tsv    = fastp_decont_report.tsv
		File?        tbprof_bam_all_depths      = collate_bam_depth.outfile
		File?        tbprof_bam_all_strains     = collate_bam_strains.outfile
		File?        tbprof_bam_all_resistances = collate_bam_resistance.outfile
		File?        tbprof_fq_all_depths       = collate_fq_depth.outfile
		File?        tbprof_fq_all_strains      = collate_fq_strains.outfile
		File?        tbprof_fq_all_resistances  = collate_fq_resistance.outfile

		# not as many stats as myco_raw for now as I'm still figuring out the best way to handle this in the multi-sample case
		# these are also not coerced, unlike the myco_raw versions, so I'm not sure if that means myco_raw and myco_sra will
		# output different things in the not defined case.
		#Float   tbd_qc_q30_in                 = fastp_decontam_check.q30_in[0]
		#Float   tbd_qc_pct_reads_NTM          = fastp_decontam_check.pct_reads_NTM[0] # TODO: how to handle for varpipe?
		#Int     tbd_qc_reads_adapter_trimmed  = fastp_decontam_check.reads_adapter_trimmed[0]
		Float?  tbd_qc_median_depth_per_tbprof = single_sample_tbprof_fq_median_depth
		Float?  tbd_qc_avg_depth_per_tbprof    = single_sample_tbprof_fq_avg_depth
		Float?  tbd_qc_pct_mapped_per_tbprof   = single_sample_tbprof_fq_pct_mapped
		Float?  tbd_qc_pct_genome_covered      = single_sample_tbprof_fq_pct_covered
		Int?    tbd_n_dr_variants              = single_sample_tbprof_fq_n_dr_variants
		Int?    tbd_n_other_variants           = single_sample_tbprof_fq_n_other_variants
		String? tbd_resistance                 = single_sample_tbprof_fq_resistance
		String? tbd_strain_per_tbprof          = single_sample_tbprof_fq_strain
	}
}

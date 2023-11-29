version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/fastp/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.11.3/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.17/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/tree_nine/0.0.11/tree_nine.wdl" as build_treesWF
import "https://raw.githubusercontent.com/aofarrel/parsevcf/1.2.0/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.2.2/tbprofiler_tasks.wdl" as profiler
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/main/thiagen_tbprofiler.wdl" as qc_fastqsWF # aka earlyQC
import "https://raw.githubusercontent.com/aofarrel/goleft-wdl/0.1.2/goleft_functions.wdl" as goleft

workflow myco {
	input {
		Array[Array[File]] paired_fastq_sets
		
		Int     covstatsQC_min_coverage        =   10
		Int     covstatsQC_min_pct_unmapped    =    2
		Boolean covstatsQC_skip_entirely       = false
		Boolean decontam_use_CDC_varpipe_ref   = true
		File?   diffQC_mask_bedfile
		Int     diffQC_max_pct_low_coverage    =    20
		Int     diffQC_this_is_low_coverage    =    10
		Int     QC_min_q30                     =    90
		Boolean clean_before_decontam          = true
		Boolean clean_after_decontam           = false
		Int     clean_average_q_score          = 29
		Boolean guardrail_mode                 = false
		Boolean soft_pct_mapped                = false
		Int     quick_tasks_disk_size          =   10 
		Boolean tbprofiler_on_bam              = false
		Int     tbprofilerQC_min_pct_mapped    =   98
		Boolean tree_decoration                = false
		File?   tree_to_decorate
	}

	parameter_meta {
		covstatsQC_min_coverage: "If covstats thinks MEAN coverage is below this, throw out this sample - not to be confused with TBProfiler MEDIAN coverage"
		covstatsQC_min_pct_unmapped: "If covstats thinks less than this percent (as int, 50 = 50%) of data does not map to H37Rv, throw out this sample"
		covstatsQC_skip_entirely: "Should we skip covstats entirely?"
		diffQC_mask_bedfile: "Bed file of regions to mask when making diff files (default: R00000039_repregions.bed)"
		diffQC_max_pct_low_coverage: "Samples who have more than this percent (as int, 50 = 50%) of positions with coverage below diffQC_this_is_low_coverage will be discarded"
		diffQC_this_is_low_coverage: "Positions with coverage below this value will be masked in diff files"
		earlyQC_min_q30_rate: "Decontaminated samples with less than this percent (as int, 50 = 50%) of reads above qual score of 30 will be discarded."
		clean_average_q_score: "Trim reads with an average quality score below this value. Independent of earlyQC_min_q30_rate. Overridden by clean_before_decontam and clean_after_decontam BOTH being false."
		quick_tasks_disk_size: "Disk size in GB to use for quick file-processing tasks; increasing this might slightly speed up file localization"
		paired_fastq_sets: "Nested array of paired fastqs, each inner array representing one samples worth of paired fastqs"
		tbprofiler_on_bam: "If true, run TBProfiler on BAMs"
		tbprofilerQC_min_pct_mapped: "If tbprofiler thinks less than this percent (as int, 50 = 50%) of data does map to H37Rv, throw out this sample"
		tree_decoration: "Should usher, taxonium, and NextStrain trees be generated?"
		tree_to_decorate: "Base tree to use if tree_decoration = true"
	}
											  
	String pass = "PASS" # used later... much later
	
	# flip some QC stuff around
	Float diffQC_max_pct_low_coverage_float = diffQC_max_pct_low_coverage / 100.0
	Int covstatsQC_max_percent_unmapped = 100 - covstatsQC_min_pct_unmapped

	scatter(paired_fastqs in paired_fastq_sets) {
		call clckwrk_combonation.clean_and_decontam_and_check as decontam_each_sample {
			input:
				docker_image = if decontam_use_CDC_varpipe_ref then "ashedpotatoes/clockwork-plus:v0.11.3.9-CDC" else "ashedpotatoes/clockwork-plus:v0.11.3.9-full",
				unsorted_sam = true,
				reads_files = paired_fastqs,
				fastp_clean_avg_qual = clean_average_q_score,
				fastp_clean_before_decontam = clean_before_decontam,
				fastp_clean_after_decontam = clean_after_decontam,
				QC_min_q30 = QC_min_q30 / 100.0,
				preliminary_min_q30 = if guardrail_mode then 0.2 else 0.0000001,
				subsample_cutoff = if guardrail_mode then 30000 else -1,
				timeout_map_reads = if guardrail_mode then 300 else 0,
				timeout_decontam = if guardrail_mode then 600 else 0
		}

		if(defined(decontam_each_sample.decontaminated_fastq_1)) {
			# This region only executes if decontaminated fastqs exist. We can use this to coerce File? into File by using
			# select_first() where the first element is the File? we know must exist, and the second element is bogus.
    		File real_decontaminated_fastq_1=select_first([decontam_each_sample.decontaminated_fastq_1, paired_fastqs[0]])
    		File real_decontaminated_fastq_2=select_first([decontam_each_sample.decontaminated_fastq_2, paired_fastqs[0]])
    		
			call qc_fastqsWF.ThiagenTBProfiler as qc_fastqs {
				input:
					fastq1 = real_decontaminated_fastq_1,
					fastq2 = real_decontaminated_fastq_2,
					soft_pct_mapped = soft_pct_mapped,
					soft_coverage = if guardrail_mode then false else true,
					minimum_coverage = if guardrail_mode then 3 else 0,
					minimum_pct_mapped = tbprofilerQC_min_pct_mapped,
					sample = decontam_each_sample.sample
			}
			# if this sample passes...
			if(qc_fastqs.status_code == pass) {
				File possibly_fastp_cleaned_fastq1=select_first([decontam_each_sample.decontaminated_fastq_1, real_decontaminated_fastq_1])
		    	File possibly_fastp_cleaned_fastq2=select_first([decontam_each_sample.decontaminated_fastq_1, real_decontaminated_fastq_2])
				call clckwrk_var_call.variant_call_one_sample_ref_included as variant_calling {
					input:
						reads_files = [possibly_fastp_cleaned_fastq1, possibly_fastp_cleaned_fastq2],
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
			
			if(covstats.percentUnmapped < covstatsQC_max_percent_unmapped) {
				if(covstats.coverage > covstatsQC_min_coverage) {
					
					# make diff files
					call diff.make_mask_and_diff as make_mask_and_diff_after_covstats {
						input:
							bam = vcfs_and_bams.left[0],
							vcf = vcfs_and_bams.right,
							min_coverage_per_site = diffQC_this_is_low_coverage,
							tbmf = diffQC_mask_bedfile,
							max_ratio_low_coverage_sites_per_sample = diffQC_max_pct_low_coverage_float
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
					min_coverage_per_site = diffQC_this_is_low_coverage,
					tbmf = diffQC_mask_bedfile,
					max_ratio_low_coverage_sites_per_sample = diffQC_max_pct_low_coverage_float
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
	Array[String] coerced_bam_strains=select_all(profile_bam.strain)
	Array[String] coerced_bam_resistances=select_all(profile_bam.resistance)
	Array[Int]    coerced_bam_depths=select_all(profile_bam.median_depth_as_int)
	
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
			Int    single_sample_tbprof_bam_depth      = coerced_bam_depths[0]
			String single_sample_tbprof_bam_resistance = coerced_bam_resistances[0]
			String single_sample_tbprof_bam_strain     = coerced_bam_strains[0]
		}
	}
  	
  	# pull TBProfiler information, if we ran TBProfiler on fastqs
  	
  	# coerce optional types into required types (doesn't crash if these are null)
	Array[String] coerced_fq_strains=select_all(qc_fastqs.strain)
	Array[String] coerced_fq_resistances=select_all(qc_fastqs.resistance)
	Array[Int]    coerced_fq_depths=select_all(qc_fastqs.median_coverage)
	
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
			Int    single_sample_tbprof_fq_depth      = coerced_fq_depths[0]
			String single_sample_tbprof_fq_resistance = coerced_fq_resistances[0]
			String single_sample_tbprof_fq_strain     = coerced_fq_strains[0]
		}
	}

	if(tree_decoration) {
		if(length(real_diffs)>0) {
			# diff files must exist if tree_decoration is true, so we can force the Array[File?]?
			# into an Array[File] with the classic "select_first() with a bogus fallback" hack
			Array[File] coerced_diffs = select_first([select_all(real_diffs), minos_vcfs])
			Array[File] coerced_reports = select_first([select_all(real_reports), minos_vcfs])
			call build_treesWF.Tree_Nine as trees {
				input:
					diffs = coerced_diffs,
					input_tree = tree_to_decorate,
					coverage_reports = coerced_reports
			}
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
		if(defined(decontam_each_sample.errorcode)) {
		
			# Sidenote: If you have a task taking in decontam_each_sample.errorcode, even after this defined check, it will fail
			# command instantiation with "Cannot interpolate Array[String?] into a command string with attribute set 
			# [PlaceholderAttributeSet(None,None,None,Some( ))]".
		
			# get the first (0th) value, eg only value since there's just one sample, and coerce it into type String
			String coerced_decontam_errorcode = select_first([decontam_each_sample.errorcode[0], "WORKFLOW_ERROR_1_REPORT_TO_DEV"])
			
			# did the decontamination step return an error?
			if(!(coerced_decontam_errorcode == pass)) {          
				String decontam_ERR = coerced_decontam_errorcode
			}
		}
		
		# handle earlyQC (if it ran at all)
		
		Array[String] earlyQC_array_coerced = select_all(qc_fastqs.status_code)
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
			if (varcall_errorcode_array[0] == "PASS") { # cannot use varcall_ERR, as it is considered optional
				#if(defined(covstats.percentUnmapped)) { # this seems to always be true, unfortunately!
					if(length(covstats.percentUnmapped) > 0) {
						# there is more than zero values in the output array, so covstats must have run
						Array[Float] percentsUnmapped = select_all(covstats.percentUnmapped)
						Float percentUnmapped = percentsUnmapped[0]
						Array[Float] meanCoverages = select_all(covstats.coverage)
						Float meanCoverage = meanCoverages[0]
						
						if(percentUnmapped > covstatsQC_max_percent_unmapped) { String too_many_unmapped = "COVSTATS_LOW_PCT_MAPPED_TO_REF" 
							if(meanCoverage < covstatsQC_min_coverage) { String double_bad = "COVSTATS_BAD_MAP_AND_COVERAGE" } 
						}
						if(meanCoverage < covstatsQC_min_coverage) { String too_low_coverage = "COVSTATS_LOW_MEAN_COVERAGE" }
					}
				#}
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
	# Failed to evaluate 'warnings' (reason 1 of 1): Evaluating flatten(flatten([[select_all(qc_fastqs.warning_codes)], [select_all(warning_decontam)]])) failed: No coercion defined from wom value(s) '[["EARLYQC_88.112_PCT_ABOVE_Q30_(MIN_0.9)", "EARLYQC_99.61_PCT_MAPPED_(MIN_99.995)"]]' of type 'Array[Array[String]]' to 'Array[String]'.
	#Array[String] warnings = flatten(flatten([[select_all(qc_fastqs.warning_codes)], [select_all(warning_decontam)]]))
	
	Map[String, String] metrics_to_values = { 
		"status": select_first([finalcode, "NA"]), 
		"reads_is_contam": select_first([decontam_each_sample.reads_is_contam[0], "NA"]),  # decontamination
		"reads_reference": select_first([decontam_each_sample.reads_reference[0], "NA"]),  # decontamination
		"reads_unmapped": select_first([decontam_each_sample.reads_unmapped[0], "NA"]),    # decontamination
		"pct_above_q30": select_first([decontam_each_sample.dcntmd_pct_above_q30[0], "NA"]),                 # fastp
		"median_coverage": select_first([qc_fastqs.median_coverage[0], "NA"]),             # thiagen!TBProfiler
		"genome_pct_coverage": select_first([qc_fastqs.pct_genome_covered[0], "NA"]),      # thiagen!TBProfiler
		"mean_coverage": select_first([meanCoverage, "NA"])                                # covstats
	}
	
	call sranwrp_processing.map_to_tsv_or_csv as qc_summary {
		input:
			the_map = metrics_to_values,
			column_names = if length(paired_fastq_sets) == 1 then [basename(paired_fastq_sets[0][0])] else ["sample"],
			outfile = if length(paired_fastq_sets) == 1 then basename(paired_fastq_sets[0][0])+"_qc" else "combined_qc_report.txt"
	}
		
	output {
		# status of sample -- only valid iff this ran on only one sample
		String status_code = select_first([finalcode, pass])
		
		# raw files
		Array[File]  bais  = final_bais
		Array[File]  bams  = final_bams
		Array[File] diffs = real_diffs
		Array[File] masks = real_masks   # bedgraph
		Array[File]  vcfs  = minos_vcfs
		
		# metadata
		Array[File?] covstats_reports          = covstats.covstatsOutfile
		Array[File?] diff_reports              = real_reports
		Array[File?] tbprof_bam_jsons          = profile_bam.tbprofiler_json
		Array[File?] tbprof_bam_summaries      = profile_bam.tbprofiler_txt
		Array[File?] tbprof_fq_jsons           = qc_fastqs.tbprofiler_json
		Array[File?] tbprof_fq_looker          = qc_fastqs.tbprofiler_looker_csv
		Array[File?] tbprof_fq_laboratorian    = qc_fastqs.tbprofiler_laboratorian_report_csv
		Array[File?] tbprof_fq_lims            = qc_fastqs.tbprofiler_lims_report_csv
		
		# these outputs only exist if there are multiple samples
		File?        tbprof_bam_all_depths      = collate_bam_depth.outfile
		File?        tbprof_bam_all_strains     = collate_bam_strains.outfile
		File?        tbprof_bam_all_resistances = collate_bam_resistance.outfile
		File?        tbprof_fq_all_depths       = collate_fq_depth.outfile
		File?        tbprof_fq_all_strains      = collate_fq_strains.outfile
		File?        tbprof_fq_all_resistances  = collate_fq_resistance.outfile
		
		# these outputs only exist if we ran on a single sample
		Int?         tbprof_bam_this_depth      = single_sample_tbprof_bam_depth
		String?      tbprof_bam_this_strain     = single_sample_tbprof_bam_strain
		String?      tbprof_bam_this_resistance = single_sample_tbprof_bam_resistance
		Int?         tbprof_fq_this_depth       = single_sample_tbprof_fq_depth
		String?      tbprof_fq_this_strain      = single_sample_tbprof_fq_strain
		String?      tbprof_fq_this_resistance  = single_sample_tbprof_fq_resistance
		
		# tree nine
		File?        tree_nwk         = trees.tree_nwk
		File?        tree_usher       = trees.tree_usher_raw
		File?        tree_taxonium    = trees.tree_taxonium
		File?        tree_nextstrain  = trees.tree_nextstrain
		Array[File]? trees_nextstrain = trees.subtrees_nextstrain
		File?        distance_matrix  = trees.distance_matrix
		
		# useful debugging/run information (only valid iff this ran on only one sample)
		File qc_csv = qc_summary.tsv_or_csv
		#Array[String] pass_or_warnings = if (length(warnings) > 0) then warnings else ["PASS"]
		String? debug_decontam_ERR  = decontam_ERR
		String? debug_earlyQC_ERR   = earlyQC_ERR
		String? debug_varcall_ERR   = varcall_ERR
		String? debug_covstats_ERR  = covstats_ERR
		String? debug_vcfdiff_ERR   = vcfdiff_ERR
		Array[String]? debug_vcfdiff_errorcode_if_covstats    = vcfdiff_errorcode_if_covstats
		Array[String]? debug_vcfdiff_errorcode_if_no_covstats = vcfdiff_errorcode_if_no_covstats
		Array[String]? debug_vcfdiff_errorcode_array          = vcfdiff_errorcode_array
		Int seconds_to_map_reads = decontam_each_sample.timer_5_mapFQ[0]
		Int seconds_to_rm_contam = decontam_each_sample.timer_7_dcnFQ[0]
		String docker_used       = decontam_each_sample.docker_used[0]
	}
}

# ** The "Scattered Optional Tasks Have Weird Outputs" (SOTHWO) Issue **
# Let's say task X and Y are mutually exclusive tasks which each output a single File named outfile.
#
# scatter {
#  if (input_variable = true) { call x }
#  if (input_variable = false) { call y }
# }
#
# If input_variable is true and this is a three way scatter, it APPEARS that, from outside the scatter...
#  * x.outfile has type Array[File?] and has three Files in it
#  * y.outfile has type Array[File?] and has at least one null (not None!) in it
#  * defined(x.outfile) is true
#  * defined(y.outfile) is true (!)
#  * length(x.outfile) is 3
#  * length(y.outfile) is 0
#  * select_all(x.outfile) creates an Array[File] with three Files
#    * You can scatter on the resulting array
#    * You can output the resulting array as a workflow-level output to a Terra data table
#  * select_all(y.outfile) creates... ???????
#    * Attempting to scatter on the resulting array will simply do nothing
#    * Idk what will happen on Terra if you have an Array[File] that is ONLY null
#  * flatten(select_all(x.outfile, y.outfile)) will result in an Array[File] that DOES have at least one null
#  * flatten(select_all(x.outfile), select_all(y.outfile)) will result in an Array[File] that DOES NOT have any nulls
#    * select_all() will also coerce x.outfile and/or y.outfile from File? to File if necessary
#  * select_first(select_all(x.outfile), select_all(y.outfile)) will result in an Array[File] that DOES NOT have any nulls
#
# Implications:
#  * nulls and Nones exist in Cromwell-WDL in spite of spec implying otherwise
#  * the mere act of scattering a task defines one array for each of its outputs, and if that task never
#    runs, the arrays remain defined and full of null(s)
#  * length() does not count null values
#  * select_all() is not "recursive" so flatten(select_all(x.outfile, y.outfile)) can have null(s)
#  * Array[File] can have nulls, so it may not be much different from Array[File?] in practice
#
# So, the correct way to gather mutually exclusive scattered optional task outputs is...
# Array[File] foo = flatten(select_all(x.outfile), select_all(y.outfile))
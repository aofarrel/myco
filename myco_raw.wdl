# myco version 6.2.4-REPRO
# This is an archived version of myco that exists solely for reproducing published results -- it is HIGHLY recommended you use a more recent version!
# The "version 1.0" string below references the WDL syntax version
version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.12.2/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.12.2/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.24/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/vcf_to_diff_wdl/0.0.3/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.2.5/tbprofiler_tasks.wdl" as profiler
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.2.5/thiagen_tbprofiler.wdl" as tbprofilerFQ_WF # fka earlyQC
import "https://raw.githubusercontent.com/aofarrel/goleft-wdl/0.1.2/goleft_functions.wdl" as goleft

workflow myco {
	input {
		
		# In the original 6.2.4 of this pipeline, some defaults varied from what you see here. This is because myco_raw was orignally
		# just the CDPH version of myco_sra, and CDPH has some different standards (some of which have since changed). Because what
		# we published was always based on on myco_sra 6.2.4's defaults, we've changed the defaults of this reproducible version of
		# myco_raw to better match myco_sra 6.2.4's defaults.
		#
		# Please be aware that myco_raw currently does not support automatic downsampling, unlike myco_sra.
		
		Array[Array[File]] paired_fastq_sets
		
		String? output_sample_name
		Boolean guardrail_mode                 = true
		
		Boolean clean_after_decontam           = false
		Int     clean_average_q_score          = 29
		Boolean clean_before_decontam          = true
		Boolean covstatsQC_skip_entirely       = true   # false in original version of myco_raw 6.2.4
		Boolean decontam_use_CDC_varpipe_ref   = false  # true in original version of myco_raw 6.2.4
		File?   mask_bedfile
		Int   QC_max_pct_low_coverage_sites    =    20
		Int     QC_max_pct_unmapped            =     2
		Int     QC_min_mean_coverage           =    10
		Int     QC_min_q30                     =    90
		Boolean QC_soft_pct_mapped             = false
		Int     QC_this_is_low_coverage        =    10
		Int     quick_tasks_disk_size          =    10 
		Boolean tbprofiler_on_bam              = true  # false in original version of myco_raw 6.2.4
	}

	parameter_meta {
		paired_fastq_sets: "Nested array of paired fastqs, each inner array representing one samples worth of paired fastqs"

		output_sample_name: "Override all sample names with this string instead"
		guardrail_mode: "Implements about a half-dozen safeguards against extremely low-quality samples running for abnormally long times"

		clean_after_decontam: "Clean fqs with fastp AFTER decontaminating (not mutually exclusive with clean_before_decontam)"
		clean_average_q_score: "Trim reads with an average quality score below this value. Independent of QC_min_q30. Overridden by clean_before_decontam and clean_after_decontam BOTH being false."
		clean_before_decontam: "Clean fqs with fastp BEFORE decontamination (not mutually exclusive with clean_after_decontam)"
		covstatsQC_skip_entirely: "Should we skip covstats entirely?"
		decontam_use_CDC_varpipe_ref: "If true, use CDC varpipe decontamination reference. If false, use CRyPTIC decontamination reference."
		mask_bedfile: "Bed file of regions to mask when making diff files (default: R00000039_repregions.bed)"
		
		QC_max_pct_low_coverage_sites: "Samples who have more than this percent (as int, 50 = 50%) of positions with coverage below QC_this_is_low_coverage will be discarded"
		QC_min_mean_coverage: "If covstats thinks MEAN coverage is below this, throw out this sample - not to be confused with TBProfiler MEDIAN coverage"
		QC_max_pct_unmapped: "If covstats thinks more than this percent of your sample (after decontam and cleaning) fails to map to H37Rv, throw out this sample."
		QC_min_q30: "Decontaminated samples with less than this percent (as int, 50 = 50%) of reads above qual score of 30 will be discarded."
		QC_soft_pct_mapped: "If true, a sample failing a percent mapped check (guardrail mode's TBProfiler check and/or covstats' check as per QC_max_pct_unmapped) will throw a non-fatal warning."
		QC_this_is_low_coverage: "Positions with coverage below this value will be masked in diff files"
		
		quick_tasks_disk_size: "Disk size in GB to use for quick file-processing tasks; increasing this might slightly speed up file localization"
		tbprofiler_on_bam: "If true, run TBProfiler on BAMs"
	}
											  
	String pass = "PASS" # used later... much later
	Float QC_max_pct_low_coverage_sites_float = QC_max_pct_low_coverage_sites / 100.0

	scatter(paired_fastqs in paired_fastq_sets) {
		call clckwrk_combonation.clean_and_decontam_and_check as decontam_each_sample {
			input:
				docker_image = if decontam_use_CDC_varpipe_ref then "ashedpotatoes/clockwork-plus:v0.11.3.11-CDC" else "ashedpotatoes/clockwork-plus:v0.11.3.11-CRyPTIC",
				unsorted_sam = true,
				force_rename_out = output_sample_name,
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
					soft_coverage = if guardrail_mode then false else true,
					minimum_coverage = if guardrail_mode then 3 else 0,
					minimum_pct_mapped = if guardrail_mode then 10 else 0, # unlike covstats, this is a MINIMUM of % MAPPED
					sample = decontam_each_sample.sample
			}
			# if this sample passes...
			if(tbprofilerFQ.status_code == pass) {
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
	Array[Int]    coerced_bam_depths=select_all(profile_bam.sample_and_median_depth)
	
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
	Array[String] coerced_fq_strains=select_all(tbprofilerFQ.sample_and_strain)
	Array[String] coerced_fq_resistances=select_all(tbprofilerFQ.sample_and_resistance)
	Array[String] coerced_fq_depths=select_all(tbprofilerFQ.sample_and_coverage)
	
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
	
	Float this_unmapped = decontam_each_sample.reads_unmapped[0]
	Float this_kept = decontam_each_sample.reads_clck_kept[0]
	Float porp_unmapped = this_unmapped/this_kept
	Float pct_unmapped_decontam = if !clean_after_decontam then (porp_unmapped * 100) else -1.0
	
	Map[String, String] metrics_to_values = { 
		"status": select_first([finalcode, "NA"]), 
		"n_reads_contam": decontam_each_sample.reads_is_contam[0],                       # decontamination
		"n_reads_decon_reference": decontam_each_sample.reads_reference[0],              # decontamination
		"n_reads_decon_unmapped": decontam_each_sample.reads_unmapped[0],                # decontamination
		"n_reads_decon_kept": decontam_each_sample.reads_clck_kept[0],                   # decontamination
		"pct_loss_decon": decontam_each_sample.pct_loss_decon[0],                        # decontamination
		"pct_loss_cleaning": decontam_each_sample.pct_loss_cleaning[0],                  # decontamination
		"pct_mapped_tbprof": select_first([tbprofilerFQ.pct_reads_mapped[0], "NA"]),     # thiagen!TBProfiler
		"pct_unmapped_covstats": select_first([percentUnmapped, "NA"]),                  # covstats 
		"pct_unmapped_decon": pct_unmapped_decontam,                                     # decontamination
		"pct_above_q30": decontam_each_sample.dcntmd_pct_above_q30[0],                   # fastp
		"median_coverage": select_first([tbprofilerFQ.median_coverage[0], "NA"]),        # thiagen!TBProfiler
		"genome_pct_coverage": select_first([tbprofilerFQ.pct_genome_covered[0], "NA"]), # thiagen!TBProfiler
		"mean_coverage": select_first([meanCoverage, "NA"])                              # covstats
	}
	String sample_name_inputs_basename = sub(sub(sub(basename(paired_fastq_sets[0][0]), ".fastq", ""), ".gz", ""), ".fq", "")
	String sample_name_maybe_varcalled = if length(final_bams) > 0 then sub(basename(final_bams[0]), "_to_H37Rv.bam", "") else sample_name_inputs_basename
	String sample_name_maybe_manually_set = if defined(output_sample_name) then select_first([output_sample_name, "fallback"]) else sample_name_maybe_varcalled
	
	call sranwrp_processing.map_to_tsv_or_csv as qc_summary {
		input:
			the_map = metrics_to_values,
			column_names = if length(paired_fastq_sets) == 1 then [sample_name_maybe_manually_set] else ["sample"],
			outfile = if length(paired_fastq_sets) == 1 then sample_name_maybe_manually_set+"_qc" else "combined_qc_report.txt"
	}
		
	output {
		# status of sample -- only valid iff this ran on only one sample
		String status_code = select_first([finalcode, pass])
		
		# debug
		Float pct_unmapped_decontamm = pct_unmapped_decontam
		Float n_reads_decon_kept = this_kept
		Float n_reads_decon_unmapped = this_unmapped
		Float? pct_mapped_tbprof = tbprofilerFQ.pct_reads_mapped[0]
		Float? pct_unmapped_covstats = percentUnmapped
		Float? pct_loss_decon = decontam_each_sample.pct_loss_decon[0]
		Float? pct_loss_cleaning = decontam_each_sample.pct_loss_cleaning[0]
		
		# raw files
		Array[File]  bais  = final_bais
		Array[File]  bams  = final_bams
		Array[File] diffs = real_diffs
		Array[File] masks = real_masks   # bedgraph
		Array[File]  vcfs  = minos_vcfs
		
		# metadata
		Array[File?] decontam_reports          = decontam_each_sample.counts_out_tsv
		Array[File?] covstats_reports          = covstats.covstatsOutfile
		Array[File?] diff_reports              = real_reports
		Array[File?] tbprof_bam_jsons          = profile_bam.tbprofiler_json
		Array[File?] tbprof_bam_summaries      = profile_bam.tbprofiler_txt
		Array[File?] tbprof_fq_jsons           = tbprofilerFQ.tbprofiler_json
		Array[File?] tbprof_fq_looker          = tbprofilerFQ.tbprofiler_looker_csv
		Array[File?] tbprof_fq_laboratorian    = tbprofilerFQ.tbprofiler_laboratorian_report_csv
		Array[File?] tbprof_fq_lims            = tbprofilerFQ.tbprofiler_lims_report_csv
		
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


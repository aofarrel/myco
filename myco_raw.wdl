version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/count-reads/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/count-reads/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/write-csv-or-tsv-task/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/tree_nine/0.0.10/tree_nine.wdl" as build_treesWF
import "https://raw.githubusercontent.com/aofarrel/parsevcf/1.2.0/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.2.2/tbprofiler_tasks.wdl" as profiler
import "https://raw.githubusercontent.com/aofarrel/TBfastProfiler/0.0.10/neoTBfastProfiler.wdl" as qc_fastqsWF # aka earlyQC
import "https://raw.githubusercontent.com/aofarrel/goleft-wdl/0.1.2/goleft_functions.wdl" as goleft

workflow myco {
	input {
		Array[Array[File]] paired_fastq_sets
		
		Int     covstatsQC_minimum_coverage    =   10
		Int     covstatsQC_max_percent_unmapped=    2
		Boolean covstatsQC_skip_entirely       = false
		Boolean decontam_use_CDC_varpipe_ref   = true
		File?   diffQC_mask_bedfile
		Int     diffQC_max_percent_low_coverage=    20
		Int     diffQC_low_coverage_cutoff     =    10
		Int     earlyQC_minimum_percent_q30    =    90
		Boolean earlyQC_skip_entirely          = false
		Boolean earlyQC_skip_QC                = false
		Boolean earlyQC_skip_trimming          = false
		Int     earlyQC_trim_qual_below        =   30
		Int     quick_tasks_disk_size          =   10 
		Int     subsample_cutoff               =   -1
		Int     subsample_seed                 = 1965
		Boolean tbprofiler_on_bam              = false
		Float   tbprofilerQC_max_pct_unmapped  =    2
		Int     timeout_decontam_part1         =    0
		Int     timeout_decontam_part2         =    0
		Int     timeout_variant_caller         =    0
		Boolean tree_decoration                = false
		File?   tree_to_decorate
		Int     variantcalling_addl_disk       =  100
		Boolean variantcalling_crash_on_error  = false
		Boolean variantcalling_crash_on_timeout= false
		Int     variantcalling_cpu             =   16
		Boolean variantcalling_debug           = false
		Int?    variantcalling_mem_height
		Int     variantcalling_memory          =   32
		Int     variantcalling_preemptibles    =    1
		Int     variantcalling_retries         =    1
		Boolean variantcalling_ssd             = true
	}

	parameter_meta {
		covstatsQC_minimum_coverage: "If covstats thinks MEAN coverage is below this, throw out this sample - not to be confused with TBProfiler MEDIAN coverage"
		covstatsQC_max_percent_unmapped: "If covstats thinks this percent (as int, 50 = 50%) of data does not map to H37Rv, throw out this sample"
		covstatsQC_skip_entirely: "Should we skip covstats entirely?"
		diffQC_mask_bedfile: "Bed file of regions to mask when making diff files (default: R00000039_repregions.bed)"
		diffQC_max_percent_low_coverage: "Samples who have more than this percent (as int, 50 = 50%) of positions with coverage below diffQC_low_coverage_cutoff will be discarded"
		diffQC_low_coverage_cutoff: "Positions with coverage below this value will be masked in diff files"
		earlyQC_minimum_percent_q30: "Decontaminated samples with less than this percent (as int, 50 = 50%) of reads above qual score of 30 will be discarded. Negated by earlyQC_skip_QC or earlyQC_skip_entirely being false."
		earlyQC_skip_entirely: "Do not run earlyQC at all - no trimming, no QC, nothing."
		earlyQC_skip_QC: "Run earlyQC (unless earlyQC_skip_entirely is true), but do not throw out samples that fail QC. Independent of earlyQC_skip_trimming."
		earlyQC_skip_trimming: "Run earlyQC (unless earlyQC_skip_entirely is true), and remove samples that fail QC (unless earlyQC_skip_QC is true), but do not use fastp's cleaned fastqs."
		earlyQC_trim_qual_below: "Trim reads with an average quality score below this value. Independent of earlyQC_minimum_percent_q30. Overridden by earlyQC_skip_trimming or earlyQC_skip_entirely being true."
		quick_tasks_disk_size: "Disk size in GB to use for quick file-processing tasks; increasing this might slightly speed up file localization"
		paired_fastq_sets: "Nested array of paired fastqs, each inner array representing one samples worth of paired fastqs"
		subsample_cutoff: "If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)"
		subsample_seed: "Seed used for subsampling with seqtk"
		tbprofiler_on_bam: "If true, run TBProfiler on BAMs"
		tbprofilerQC_max_pct_unmapped: "If tbprofiler thinks this percent (as float, 50 = 50%) of data does not map to H37Rv, throw out this sample"
		timeout_decontam_part1: "Discard any sample that is still running in clockwork map_reads after this many minutes (set to 0 to never timeout)"
		timeout_decontam_part2: "Discard any sample that is still running in clockwork rm_contam after this many minutes (set to 0 to never timeout)"
		timeout_variant_caller: "Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to 0 to never timeout)"
		tree_decoration: "Should usher, taxonium, and NextStrain trees be generated?"
		tree_to_decorate: "Base tree to use if tree_decoration = true"
	}
											  
	String pass = "PASS" # used later... much later
	
	# convert percent integers to floats (except covstatsQC_max_percent_unmapped)
	Float diffQC_max_percent_low_coverage_float = diffQC_max_percent_low_coverage / 100.0
	Float earlyQC_minimum_percent_q30_float = earlyQC_minimum_percent_q30 / 100.0
	Float tbprofilerQC_max_percent_unmapped_float = tbprofilerQC_max_pct_unmapped / 100.0

	scatter(paired_fastqs in paired_fastq_sets) {
		call clckwrk_combonation.combined_decontamination_single_ref_included as decontam_each_sample {
			input:
				docker_image = if decontam_use_CDC_varpipe_ref then "ashedpotatoes/clockwork-plus:v0.11.3.7-CDC" else "ashedpotatoes/clockwork-plus:v0.11.3.2-full",
				unsorted_sam = true,
				reads_files = paired_fastqs,
				subsample_cutoff = subsample_cutoff,
				subsample_seed = subsample_seed,
				timeout_map_reads = timeout_decontam_part1,
				timeout_decontam = timeout_decontam_part2
		}

		if(defined(decontam_each_sample.decontaminated_fastq_1)) {
			# This region only executes if decontaminated fastqs exist. We can use this to coerce File? into File by using
			# select_first() where the first element is the File? we know must exist, and the second element is bogus.
    		File real_decontaminated_fastq_1=select_first([decontam_each_sample.decontaminated_fastq_1, paired_fastqs[0]])
    		File real_decontaminated_fastq_2=select_first([decontam_each_sample.decontaminated_fastq_2, paired_fastqs[0]])
    		
    		if((decontam_each_sample.reads_kept < 5000)) {
		        String warning_decontam = "DECONTAMINATION_ONLY" + decontam_each_sample.reads_kept + "_READS_REMAINING_(MIN_" + 5000 + ")" #!StringCoercion
		    }

			if(!earlyQC_skip_entirely) {
				call qc_fastqsWF.TBfastProfiler as qc_fastqs {
					input:
						fastq1 = real_decontaminated_fastq_1,
						fastq2 = real_decontaminated_fastq_2,
						q30_cutoff = earlyQC_minimum_percent_q30_float,
						average_qual = earlyQC_trim_qual_below,
						use_fastps_cleaned_fastqs = !(earlyQC_skip_trimming),
						soft_all_qc = earlyQC_skip_QC,
						pct_mapped_cutoff = tbprofilerQC_max_percent_unmapped_float
				}
				# if this sample passes fastp, or if earlyQC_skip_QC is true...
				if(qc_fastqs.status_code == pass) {
					File possibly_fastp_cleaned_fastq1=select_first([qc_fastqs.cleaned_fastq1, real_decontaminated_fastq_1])
			    	File possibly_fastp_cleaned_fastq2=select_first([qc_fastqs.cleaned_fastq2, real_decontaminated_fastq_2])
					call clckwrk_var_call.variant_call_one_sample_ref_included as variant_call_after_earlyQC {
						input:
							reads_files = [possibly_fastp_cleaned_fastq1, possibly_fastp_cleaned_fastq2],
							addldisk = variantcalling_addl_disk,
							cpu = variantcalling_cpu,
							crash_on_error = variantcalling_crash_on_error,
							crash_on_timeout = variantcalling_crash_on_timeout,
							debug = variantcalling_debug,
							mem_height = variantcalling_mem_height,
							memory = variantcalling_memory,
							preempt = variantcalling_preemptibles,
							retries = variantcalling_retries,
							ssd = variantcalling_ssd,
							tarball_bams_and_bais = false,
							timeout = timeout_variant_caller
					}
				}
			}
			
			# if we ARE skipping early QC (but the samples did decontaminate without erroring/timing out)
			if(earlyQC_skip_entirely) {
				call clckwrk_var_call.variant_call_one_sample_ref_included as variant_call_without_earlyQC {
					input:
						reads_files = [real_decontaminated_fastq_1, real_decontaminated_fastq_2],
						addldisk = variantcalling_addl_disk,
						cpu = variantcalling_cpu,
						crash_on_error = variantcalling_crash_on_error,
						crash_on_timeout = variantcalling_crash_on_timeout,
						debug = variantcalling_debug,
						mem_height = variantcalling_mem_height,
						memory = variantcalling_memory,
						preempt = variantcalling_preemptibles,
						retries = variantcalling_retries,
						ssd = variantcalling_ssd,
						tarball_bams_and_bais = false,
						timeout = timeout_variant_caller
				}
			}
		}

	}

	# do some wizardry to deal with optionals
	#
	# In order to account for different use cases, this workflow has two versions of the variant caller. They are mutually
	# exclusive, eg, only one can ever be called by a given sample. In fact, there are cases where NONE of them get called.
	# Additionally, each version of the variant caller technically gives optional output. This is to prevent the entire
	# pipeline from crashing if a single garbage sample starts the variant calling task but cannot make a VCF.
	#
	# Unfortunately, this means I have created multiple optional outputs from optional tasks. Mutually exclusive outputs from
	# mutually exclusive tasks, yes, but WDL doesn't quite understand mutual exclusivity. WDL and/or Cromwell (it's hard to 
	# know if the issue is the langauge or its implementation since the spec is not very specific) also gets finicky when we try
	# to do certain things with optional arrays. So, we want to turn those optionals into not-optionals ASAP. That's what this
	# block of variable definitions does. WDL does not allow you to overwrite variables, so we need to declare a ton of variables
	# with unique names. Yes, this could in theory be done more "cleanly" by calling a task, but why spin up a Docker image to do
	# microseconds worth of work?
	#
	# Interestingly, even if the variant caller didn't run, this code block does not cause a crash. Somehow, we can create 
	# non-optional variables that have no content. This little mystery is mostly good because it means the 
	# all-samples-get-dropped-before-variant-calling-because-they-suck does not break this section... but it also means
	# we cannot use anything here to test if any version of the variant caller actually ran at all!
	# 
	Array[File] minos_vcfs_if_earlyQC = select_all(variant_call_after_earlyQC.adjudicated_vcf)
	Array[File] minos_vcfs_if_no_earlyQC = select_all(variant_call_without_earlyQC.adjudicated_vcf)
	
	Array[File] bams_if_earlyQC = select_all(variant_call_after_earlyQC.bam)
	Array[File] bams_if_no_earlyQC = select_all(variant_call_without_earlyQC.bam)
	
	Array[File] bais_if_earlyQC = select_all(variant_call_after_earlyQC.bai)
	Array[File] bais_if_no_earlyQC = select_all(variant_call_without_earlyQC.bai)
	
	Array[File] minos_vcfs = flatten([minos_vcfs_if_earlyQC, minos_vcfs_if_no_earlyQC])
	Array[File] final_bams = flatten([bams_if_earlyQC, bams_if_no_earlyQC])
	Array[File] final_bais = flatten([bais_if_earlyQC, bais_if_no_earlyQC])
	
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
				if(covstats.coverage > covstatsQC_minimum_coverage) {
					
					# make diff files
					call diff.make_mask_and_diff as make_mask_and_diff_after_covstats {
						input:
							bam = vcfs_and_bams.left[0],
							vcf = vcfs_and_bams.right,
							min_coverage_per_site = diffQC_low_coverage_cutoff,
							tbmf = diffQC_mask_bedfile,
							max_ratio_low_coverage_sites_per_sample = diffQC_max_percent_low_coverage_float
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
					min_coverage_per_site = diffQC_low_coverage_cutoff,
					tbmf = diffQC_mask_bedfile,
					max_ratio_low_coverage_sites_per_sample = diffQC_max_percent_low_coverage_float
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

		# Unfortunately, to check if the variant caller ran, we have to check both versions of the variant caller.
		#
		# But there seems to be a bug in Cromwell that causes it to incorrectly define variables even when a task hasn't run.
		# For example, if(defined(variant_call_after_earlyQC.errorcode)) returns true even if NO variant callers
		# ran at all. And if you put a select_first() in that if block, it will fall back to the next variable, indicating
		# that select_first() knows it's undefined but defined() does not. I have not replicated this behavior with an optional
		# non-scattered task, so I think this is specific to scattered tasks, or I am a fool and something is wrong with my 
		# defined() checks. In any case, here's the Cromwell ticket: https://github.com/broadinstitute/cromwell/issues/7201
		#
		# I've decided to instead use a variation of how we coerce the VCF outputs into required types - relying entirely on
		# select_first() and select_all() instead of the seemingly buggy defined(). It's less clean, but it seems to be necessary.
		# The WDL 1.0 spec does not say what happens if you give select_all() an array that only has optional values, but
		# the WDL 1.1 spec says you get an empty array. Thankfully, Cromwell handles 1.0-select_all() like the 1.1 spec.
		Array[String] errorcode_if_earlyQC = select_all(variant_call_after_earlyQC.errorcode)
		Array[String] errorcode_if_no_earlyQC = select_all(variant_call_without_earlyQC.errorcode)
		
		# if the variant caller did not run, the fallback pass will be selected, even though the sample shouldn't be considered a pass, so
		# the final-final-final error code needs to have decontam's error come before the variant caller error.
		Array[String] varcall_errorcode_array = flatten([errorcode_if_earlyQC, errorcode_if_no_earlyQC, ["PASS"]])
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
							if(meanCoverage < covstatsQC_minimum_coverage) { String double_bad = "COVSTATS_BAD_MAP_AND_COVERAGE" } 
						}
						if(meanCoverage < covstatsQC_minimum_coverage) { String too_low_coverage = "COVSTATS_LOW_MEAN_COVERAGE" }
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
	Array[String] warnings = flatten(flatten([[select_all(qc_fastqs.warning_codes)], [select_all(warning_decontam)]]))
	
	call sranwrp_processing.write_csv as qc_summary {
		input:
			# TODO: warnings has commas, maybe don't include that? or see if TSV will work for CDPH
			headings = [
				"status",
				"reads_is_contam", 
				"reads_reference",
				"reads_unmapped", 
				"median_coverage", 
				"mean_coverage" ],
			stuff_to_write = [[
				select_first([finalcode, "NA"]),
				select_first([decontam_each_sample.reads_is_contam[0], "NA"]), 
				select_first([decontam_each_sample.reads_reference[0], "NA"]), 
				select_first([decontam_each_sample.reads_unmapped[0], "NA"]),
				select_first([qc_fastqs.median_coverage[0], single_sample_tbprof_fq_depth, "NA"]),
				select_first([meanCoverage, "NA"]) ]]
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
		Array[File?] fastp_reports             = qc_fastqs.fastp_txt
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
		
		# useful debugging/run information (only valid iff this ran on only one sample)
		File qc_csv = qc_summary.finalOut
		Array[String] pass_or_warnings = if (length(warnings) > 0) then warnings else ["PASS"]
		String? debug_decontam_ERR  = decontam_ERR
		String? debug_earlyQC_ERR   = earlyQC_ERR
		String? debug_varcall_ERR   = varcall_ERR
		String? debug_covstats_ERR  = covstats_ERR
		String? debug_vcfdiff_ERR   = vcfdiff_ERR
		Array[String]? debug_vcfdiff_errorcode_if_covstats    = vcfdiff_errorcode_if_covstats
		Array[String]? debug_vcfdiff_errorcode_if_no_covstats = vcfdiff_errorcode_if_no_covstats
		Array[String]? debug_vcfdiff_errorcode_array          = vcfdiff_errorcode_array
		Int seconds_to_map_reads = decontam_each_sample.timer_map_reads[0]
		Int seconds_to_rm_contam = decontam_each_sample.timer_rm_contam[0]
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
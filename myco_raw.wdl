version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/error-codes/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/error-codes/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.12/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/tree_nine/0.0.10/tree_nine.wdl" as build_treesWF
import "https://raw.githubusercontent.com/aofarrel/parsevcf/1.1.9/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.2.2/tbprofiler_tasks.wdl" as profiler
import "https://raw.githubusercontent.com/aofarrel/TBfastProfiler/record-if-sample-fails-qc/TBfastProfiler.wdl" as qc_fastqsWF # aka earlyQC
import "https://raw.githubusercontent.com/aofarrel/goleft-wdl/0.1.2/goleft_functions.wdl" as goleft


workflow myco {
	input {
		Array[Array[File]] paired_fastq_sets
		
		# TODO: REPLACE WITH BETTER DEFAULTS
		Float   covstats_qc_cutoff_coverages    =   10.00
		Float   covstats_qc_cutoff_unmapped     =    2.00
		Boolean covstats_qc_skip_entirely       = false
		
		Boolean diff_force                      = false
		File?   diff_mask_these_regions
		Int     diff_min_coverage_per_site      =   10
		Boolean early_qc_apply_cutoffs          = false
		Float   early_qc_cutoff_q30             =    0.90
		Boolean early_qc_skip_entirely          = false
		Int     quick_tasks_disk_size           =   10 
		Int     subsample_cutoff                =   -1
		Int     subsample_seed                  = 1965
		Boolean tbprofiler_on_bam               = false
		Int     timeout_decontam_part1          =    0
		Int     timeout_decontam_part2          =    0
		Int     timeout_variant_caller          =    0
		Boolean tree_decoration                 = false
		Float   tree_max_low_coverage_sites     =    0.05
		File?   tree_to_decorate
		Int     variantcalling_addl_disk        =  100
		Boolean variantcalling_crash_on_error   = false
		Boolean variantcalling_crash_on_timeout = false
		Int     variantcalling_cpu              =   16
		Boolean variantcalling_debug            = false
		Int?    variantcalling_mem_height
		Int     variantcalling_memory           =   32
		Int     variantcalling_preemptibles     =    1
		Int     variantcalling_retries          =    1
		Boolean variantcalling_ssd              = true
	}

	parameter_meta {
		covstats_qc_cutoff_coverages: "If covstats thinks coverage is below this, throw out this sample"
		covstats_qc_cutoff_unmapped: "If covstats thinks this porportion (as float, 50 = 50%) of data does not map to H37Rv, throw out this sample"
		covstats_qc_skip_entirely: "Should we skip covstats entirely?"
		diff_force: "If true and if decorate_tree is false, generate diff files. (Diff files will always be created if decorate_tree is true.)"
		diff_mask_these_regions: "Bed file of regions to mask when making diff files"
		diff_min_coverage_per_site: "Positions with coverage below this value will be masked in diff files"
		early_qc_apply_cutoffs: "If true, run fastp + TBProfiler on decontaminated fastqs and apply cutoffs to determine which samples should be thrown out."
		early_qc_cutoff_q30: "Decontaminated samples with less than this porportion (as float, 0.5 = 50%) of reads above qual score of 30 will be discarded iff early_qc_apply_cutoffs is also true."
		early_qc_skip_entirely: "Do not run early QC (fastp + fastq-TBProfiler) at all. Does not affect whether or not TBProfiler is later run on bams. Overrides early_qc_apply_cutoffs."
		quick_tasks_disk_size: "Disk size in GB to use for quick file-processing tasks; increasing this might slightly speed up file localization"
		paired_fastq_sets: "Nested array of paired fastqs, each inner array representing one samples worth of paired fastqs"
		ref_genome_for_tree_building: "Ref genome for building trees -- must have ONLY `>NC_000962.3` on its first line"
		subsample_cutoff: "If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)"
		subsample_seed: "Seed used for subsampling with seqtk"
		tbprofiler_on_bam: "If true, run TBProfiler on BAMs"
		timeout_decontam_part1: "Discard any sample that is still running in clockwork map_reads after this many minutes (set to 0 to never timeout)"
		timeout_decontam_part2: "Discard any sample that is still running in clockwork rm_contam after this many minutes (set to 0 to never timeout)"
		timeout_variant_caller: "Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to 0 to never timeout)"
		tree_decoration: "Should usher, taxonium, and NextStrain trees be generated?"
		tree_max_low_coverage_sites: "If a diff file has higher than this porportion (0.5 = 50%) bad data, do not include it in the tree"
		tree_to_decorate: "Base tree to use if tree_decoration = true"
	}

	# WDL doesn't understand mutual exclusivity, so we have to get a little creative on 
	# our determination of whether or not we want to create diff files.
	if(tree_decoration)  {  Boolean create_diff_files_   = true  }
	if(!tree_decoration) {
		if(!tree_decoration){  Boolean create_diff_files__  = false }
		if(tree_decoration) {  Boolean create_diff_files___ = true  }
	}
	Boolean create_diff_files = select_first([create_diff_files_,
											  create_diff_files__, 
											  create_diff_files___])
											  
	String pass = "PASS" # used later... much later

	scatter(paired_fastqs in paired_fastq_sets) {
		call clckwrk_combonation.combined_decontamination_single_ref_included as decontam_each_sample {
			input:
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

			if(!early_qc_skip_entirely) {
				call qc_fastqsWF.TBfastProfiler as qc_fastqs {
					input:
						fastq1 = real_decontaminated_fastq_1,
						fastq2 = real_decontaminated_fastq_2,
						q30_cutoff = early_qc_cutoff_q30
				}
				
				# if we are filtering out samples via earlyQC...
				if(early_qc_apply_cutoffs) {
				
					# and this sample passes...
					if(qc_fastqs.pass_or_errorcode == pass) {
						File possibly_fastp_cleaned_fastq1_passed=select_first([qc_fastqs.cleaned_fastq1, real_decontaminated_fastq_1])
				    	File possibly_fastp_cleaned_fastq2_passed=select_first([qc_fastqs.cleaned_fastq2, real_decontaminated_fastq_2])
						call clckwrk_var_call.variant_call_one_sample_ref_included as variant_call_after_earlyQC_filtering {
							input:
								reads_files = [possibly_fastp_cleaned_fastq1_passed, possibly_fastp_cleaned_fastq2_passed],
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
				
				# if we are not filtering out samples via the early qc step (but ran earlyQC anyway)...
				if(!early_qc_apply_cutoffs) {
					File possibly_fastp_cleaned_fastq1=select_first([qc_fastqs.cleaned_fastq1, real_decontaminated_fastq_1])
			    	File possibly_fastp_cleaned_fastq2=select_first([qc_fastqs.cleaned_fastq2, real_decontaminated_fastq_2])
					call clckwrk_var_call.variant_call_one_sample_ref_included as variant_call_after_earlyQC_but_not_filtering_samples {
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
			if(early_qc_skip_entirely) {
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
	# In order to account for different use cases, this workflow has three versions of the variant caller. They are mutually
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
	Array[File] minos_vcfs_if_earlyQC_filtered = select_all(variant_call_after_earlyQC_filtering.adjudicated_vcf)
	Array[File] minos_vcfs_if_earlyQC_but_not_filtering = select_all(variant_call_after_earlyQC_but_not_filtering_samples.adjudicated_vcf)
	Array[File] minos_vcfs_if_no_earlyQC = select_all(variant_call_without_earlyQC.adjudicated_vcf)
	
	Array[File] bams_if_earlyQC_filtered = select_all(variant_call_after_earlyQC_filtering.bam)
	Array[File] bams_if_earlyQC_but_not_filtering = select_all(variant_call_after_earlyQC_but_not_filtering_samples.bam)
	Array[File] bams_if_no_earlyQC = select_all(variant_call_without_earlyQC.bam)
	
	Array[File] bais_if_earlyQC_filtered = select_all(variant_call_after_earlyQC_filtering.bai)
	Array[File] bais_if_earlyQC_but_not_filtering = select_all(variant_call_after_earlyQC_but_not_filtering_samples.bai)
	Array[File] bais_if_no_earlyQC = select_all(variant_call_without_earlyQC.bai)
	
	Array[File] minos_vcfs = flatten([minos_vcfs_if_earlyQC_filtered, minos_vcfs_if_earlyQC_but_not_filtering, minos_vcfs_if_no_earlyQC])
	Array[File] final_bams = flatten([bams_if_earlyQC_filtered, bams_if_earlyQC_but_not_filtering, bams_if_no_earlyQC])
	Array[File] final_bais = flatten([bais_if_earlyQC_filtered, bais_if_earlyQC_but_not_filtering, bais_if_no_earlyQC])
	
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
		if(!covstats_qc_skip_entirely) {
	
			# covstats to check coverage and percent mapped to reference
			call goleft.covstats as covstats {
				input:
					inputBamOrCram = vcfs_and_bams.left[0],
					allInputIndexes = [vcfs_and_bams.left[1]]
			}
			
			if(covstats.percentUnmapped > covstats_qc_cutoff_unmapped) {
				if(covstats.coverage > covstats_qc_cutoff_coverages) {
					
					# make diff files
					call diff.make_mask_and_diff as make_mask_and_diff_after_covstats {
						input:
							bam = vcfs_and_bams.left[0],
							vcf = vcfs_and_bams.right,
							min_coverage_per_site = diff_min_coverage_per_site,
							tbmf = diff_mask_these_regions,
							diffs = create_diff_files
					}
				}
			}
		}
		
		if(covstats_qc_skip_entirely) {
		
			# make diff files
			call diff.make_mask_and_diff as make_mask_and_diff_no_covstats {
				input:
					bam = vcfs_and_bams.left[0],
					vcf = vcfs_and_bams.right,
					min_coverage_per_site = diff_min_coverage_per_site,
					tbmf = diff_mask_these_regions,
					diffs = create_diff_files
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
	
	Array[File?] real_diffs = select_first([make_mask_and_diff_after_covstats.diff, make_mask_and_diff_no_covstats.diff])
	Array[File?] real_reports = select_first([make_mask_and_diff_after_covstats.report, make_mask_and_diff_no_covstats.report])
	Array[File?] real_masks = select_first([make_mask_and_diff_after_covstats.mask_file, make_mask_and_diff_no_covstats.mask_file])

	# pull TBProfiler information, if we ran TBProfiler on bams
	# For reasons I don't understand, if no variant caller runs (ergo there's no bam and profile_bam also does not run),
	# the following code block will still execute. coerced_bam_strains etc will be arrays with length 0. So this check
	# isn't sufficient on its own, but as always, I have a workaround for this!
	if(defined(profile_bam.strain)) {
	
		# coerce optional types into required types
		Array[String] coerced_bam_strains=select_all(profile_bam.strain)
		Array[String] coerced_bam_resistances=select_all(profile_bam.resistance)
		Array[String] coerced_bam_depths=select_all(profile_bam.median_depth)
		
		# workaround for "profile_bam.strain exists but profile_bam didn't run" bug
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
				Int    single_sample_tbprof_bam_depth      = read_int(coerced_bam_depths[0])
				String single_sample_tbprof_bam_resistance = coerced_bam_resistances[0]
				String single_sample_tbprof_bam_strain     = coerced_bam_strains[0]
			}
		}
  	}
  	
  	# pull TBProfiler information, if we ran TBProfiler on fastqs
  	if(defined(qc_fastqs.samp_strain)) {
		Array[String] coerced_fq_strains=select_all(qc_fastqs.samp_strain)
		Array[String] coerced_fq_resistance=select_all(qc_fastqs.samp_resistance)

		call sranwrp_processing.cat_strings as collate_fq_strains {
			input:
				strings = coerced_fq_strains,
				out = "strain_reports.txt",
				disk_size = quick_tasks_disk_size
		}
		
		call sranwrp_processing.cat_strings as collate_fq_resistance {
			input:
				strings = coerced_fq_resistance,
				out = "resistance_reports.txt",
				disk_size = quick_tasks_disk_size
		}
  	}

	if(tree_decoration) {
		# diff files must exist if tree_decoration is true, so we can force the Array[File?]?
		# into an Array[File] with the classic "select_first() with a bogus fallback" hack
		Array[File] coerced_diffs = select_first([select_all(real_diffs), minos_vcfs])
		Array[File] coerced_reports = select_first([select_all(real_reports), minos_vcfs])
		call build_treesWF.Tree_Nine as trees {
			input:
				diffs = coerced_diffs,
				input_tree = tree_to_decorate,
				coverage_reports = coerced_reports,
				max_low_coverage_sites = tree_max_low_coverage_sites
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
	
		# did the decontamination step actually run?
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
		
		# We already set earlyqc_definitely_errored earlier, but since this workflow can run on multiple samples it's going to have
		# type Array[String?], Array[String?]?, or Array[String]? (should probably be Array[String?]? but optionals are so messy it
		# can be hard to predict what the interpreter thinks).
		
		if (defined(qc_fastqs.pass_or_errorcode)) {
			String coerced_earlyqc_errorcode = select_first([qc_fastqs.pass_or_errorcode[0], "WORKFLOW_ERROR_999_REPORT_TO_DEV"])
			
			# did the sample pass earlyqc?
			if(!(coerced_decontam_errorcode == pass)) {          
				String earlyqc_ERR = coerced_earlyqc_errorcode
			}
		}

		# Unfortunately, to check if the variant caller ran, we have to check all three versions of the variant caller.	We also
		# have to use the bogus select_first() fallback to define the new strings, because WDL doesn't understand that if you 
		# are in a block that only executes if X is defined, then X must be defined. (This isn't a Cromwell thing; miniwdl 
		# also catches this.)
		#
		# But what miniwdl (although this may be a Cromwell-runtime thing, I'm not sure) and womtool don't find is this: we gotta use
		# the bogus select_first() workaround not only to define our new strings, but also the statement immediately prior
		# that checks if the errorcode is "PASS", or else at runtime Cromwell will return "Failed to evaluate 'if_condition'", b/c
		# "Evaluating !((variant_call_without_earlyQC.errorcode[0] == pass)) failed: Sorry! Operation == is not supported on
		# empty optional values. You might resolve this using select_first([optional, default]) to guarantee that you have a
		# filled value."" I'm not sure why this can't be detected before runtime.
		#
		
		# TODO: this keeps falling back to error code 2. I'm hoping it's a call cache bug?
		# did the "if earlyQC filtered" variant caller run?
		if(defined(variant_call_after_earlyQC_filtering.errorcode)) {

			# I think what's happening is that the variable is incorrectly being considered defined when it shouldn't be, since I get 
			# error two when running with skipping early qc. so varcall_ERR selects varcall_error_if_earlyQC_filtered's fallback, even though
			# varcall_error_if_no_earlyQC is valid...
			
			if(length(variant_call_after_earlyQC_filtering.errorcode) > 0) {
				# get the first (0th) value and coerce it into type String
				String coerced_vc_filtered_errorcode = select_first([variant_call_after_earlyQC_filtering.errorcode[0], "WORKFLOW_ERROR_2_REPORT_TO_DEV"])
				call echo {input: integer=length(variant_call_after_earlyQC_filtering.errorcode), string=variant_call_after_earlyQC_filtering.errorcode[0]}
				call echo_array {input: array=[coerced_vc_filtered_errorcode]}
			}
			if(!(length(variant_call_after_earlyQC_filtering.errorcode) > 0)) {
				# get the first (0th) value and coerce it into type String
				String coerced_vc_filtered_errorcode_alt = select_first([variant_call_after_earlyQC_filtering.errorcode[0], "WORKFLOW_ERROR_3_REPORT_TO_DEV"])
			}
			
			# did the "if earlyQC filtered" variant caller return an error?
			String foo = select_first([coerced_vc_filtered_errorcode_alt, coerced_vc_filtered_errorcode, "ARGH"])
			if(!(foo == pass)) {
				String varcall_error_if_earlyQC_filtered = select_first([coerced_vc_filtered_errorcode_alt, coerced_vc_filtered_errorcode, "AAAA"])

			}
		}
		# did the "if earlyQC but not filtered" variant caller run?
		if(defined(variant_call_after_earlyQC_but_not_filtering_samples.errorcode)) {
		
			# get the first (0th) value and coerce it into type String
			String coerced_vc_notfiltered_errorcode = select_first([variant_call_after_earlyQC_but_not_filtering_samples.errorcode[0], "WORKFLOW_ERROR_4_REPORT_TO_DEV"])
			
			# did the "if earlyQC but not filtered" variant caller return an error?		
			# TODO: do we need this select_first()?
			if(!(select_first([coerced_vc_notfiltered_errorcode, "silly bogus fallback"]) == pass)) {
				String varcall_error_if_earlyQC_but_not_filtering = select_first([coerced_vc_notfiltered_errorcode, "WORKFLOW_ERROR_5_REPORT_TO_DEV"])
			}
		}
		# did the "no earlyQC" variant caller run?
		if(defined(variant_call_without_earlyQC.errorcode)) {
		
			# get the first (0th) value and coerce it into type String
			String coerced_vc_noearlyqc_errorcode = select_first([variant_call_without_earlyQC.errorcode[0], "WORKFLOW_ERROR_6_REPORT_TO_DEV"])
		
			# did the "no earlyQC" variant caller return an error?
			# TODO: do we need select_first()?
			if(!(select_first([coerced_vc_noearlyqc_errorcode, "silly bogus fallback"]) == pass)) {
				String varcall_error_if_no_earlyQC = select_first([coerced_vc_noearlyqc_errorcode, "WORKFLOW_ERROR_7_REPORT_TO_DEV"])
			}
		}
		
		# if the variant caller did not run, the fallback is pass, even though the sample shouldn't be considered a pass, so
		# the final-final-final error code needs to have decontam's error come before the variant caller error.
		String varcall_ERR = select_first([varcall_error_if_earlyQC_filtered, varcall_error_if_earlyQC_but_not_filtering, varcall_error_if_no_earlyQC, pass])
	
		# final-final-final error code
		String finalcode = select_first([decontam_ERR, earlyqc_ERR, varcall_ERR, pass])

	}
		
	output {
		# raw files
		Array[File]  bais = final_bais
		Array[File]  bams  = final_bams
		Array[File?] diffs = real_diffs
		Array[File?] masks = real_masks   # bedgraph
		Array[File]  vcfs  = minos_vcfs
		
		# metadata
		Array[File?]  covstats_reports       = covstats.covstatsOutfile
		Array[File?]  diff_reports           = real_reports
		Array[File?]  fastp_reports          = qc_fastqs.fastp_txt
		Array[File?]  tbprof_bam_jsons       = profile_bam.tbprofiler_json
		Array[File?]  tbprof_bam_summaries   = profile_bam.tbprofiler_txt
		Array[File?]  tbprof_fq_jsons        = qc_fastqs.tbprofiler_json
		File?         tbprof_fq_strains      = collate_fq_strains.outfile
		Array[File?]  tbprof_fq_summaries    = qc_fastqs.tbprofiler_txt
		File?         tbprof_fq_resistances  = collate_fq_resistance.outfile
		
		# these outputs only exist if there are multiple samples
		File?         tbprof_bam_depths      = collate_bam_depth.outfile
		File?         tbprof_bam_strains     = collate_bam_strains.outfile
		File?         tbprof_bam_resistances = collate_bam_resistance.outfile
		
		# these outputs only exist if we ran on a single sample
		Int?            tbprof_bam_depth      = single_sample_tbprof_bam_depth
		String?         tbprof_bam_strain     = single_sample_tbprof_bam_strain
		String?         tbprof_bam_resistance = single_sample_tbprof_bam_resistance
		
		# status of sample, only valid iff this ran on only one sample
		String error_code = select_first([finalcode, pass])
		
		# tree nine
		File?        tree_nwk         = trees.tree_nwk
		File?        tree_usher       = trees.tree_usher_raw
		File?        tree_taxonium    = trees.tree_taxonium
		File?        tree_nextstrain  = trees.tree_nextstrain
		Array[File]? trees_nextstrain = trees.subtrees_nextstrain
	}
}

task echo {
	input {
		Int? integer
		String? string
	}
	
	command  <<<
	echo "~{integer}"
	echo "~{string}"
	>>>
	
	runtime {
		cpu: 2
		docker: "ashedpotatoes/iqbal-unofficial-clockwork-mirror:v0.11.3"
		disks: "local-disk " + 10 + " HDD"
		memory: "4 GB"
		preemptible: "1"
	}
}

task echo_array {
	input {
		Array[String?] array
	}
	
	command  <<<
	echo "~{sep=' ' array}"
	>>>
	
	runtime {
		cpu: 2
		docker: "ashedpotatoes/iqbal-unofficial-clockwork-mirror:v0.11.3"
		disks: "local-disk " + 10 + " HDD"
		memory: "4 GB"
		preemptible: "1"
	}
}
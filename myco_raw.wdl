version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.16.9/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.16.9/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.24/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/vcf_to_diff_wdl/0.0.5/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.3.2/tbprofiler_tasks.wdl" as profiler
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.3.2/theiagen_tbprofiler.wdl" as tbprofilerFQ_WF # fka earlyQC
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
		String output_sample_name

		File?   call_as_reference_bedfile           # default: R00000039_repregions.bed (exists in the Docker image)
		String? comment
		Int     fastp_avg_qual              = 29
		Boolean just_like_2024              = false
		Boolean guardrail_mode              = true
		Boolean low_resource_mode           = false
		Int     sample_max_pct_masked       = 20
		Int     sample_min_pct_mapped       = 90
		Int     sample_min_avg_depth        = 30
		Int     sample_min_q30              = 80
		Int     site_min_depth              = 10
		Boolean skip_covstats               = true
		Int     subsample_cutoff            = -1     # note inconsistency with myco_sra and how guardrail_mode affects this
	}

	parameter_meta {
		# parameter_meta doesn't support multi-line strings, so I'll include additional information in comments to avoid making this file super wide
		#
		# There's basically three levels of QC going on here:
		# * pass/fail an entire sample (note that "fail" does NOT mean "pipeline will return 1")
		# * pass/fail an entire read/pair, where failing reads are trimmed or discarded as per fastp standards
		# * pass/fail a specific site, where failing sites are masked in the final diff file

		call_as_reference_bedfile: "Bed file of regions to mask to reference when making diff files (default: R00000039_repregions.bed)"
		# Default: https://github.com/iqbal-lab-org/cryptic_tb_callable_mask/blob/44f884558bea4ee092ce7c5c878561200fcee92f/R00000039_repregions.bed
		# Note that masking *to reference* is not the same as low-coverage masking. Reference (and masked-to reference) positions are not mentioned
		# in diff files at all, but when we low-coverage mask we explictly say -. Reference and - may be treated differently by matUtils/UShER.

		comment: "Information about this run; this is echoed as a workflow level output"
		# This is useful for tracking provenance/versioning for Terra data tables
		
		fastp_avg_qual: "Discard read pairs with an average quality score below this value via fastp --average_qual. Not the same as sample_min_q30."
		# Fastp defines this as:
		# "if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])"
		# What you set this value to you might affect the results of sample_min_q30, but keep in mind that is a whole-sample filter.

		#decontam_use_CDC_varpipe_ref: "If true, use CDC varpipe decontamination reference. If false, use CRyPTIC decontamination reference."
		# CDC uses their own version of clockwork's decontamination reference, which I call "CDC varpipe" since I pulled it from the varpipe repo.
		# This is currently a null op since it'd require I maintain double the number of Docker images, and it doesn't delineate between human vs
		# NTM vs other forms of contamination (which the task currently requires for some outputs). If there is a demand for CDC varpipe I can
		# make this an option again, but I nevertheless gently recommend against using it due to unclear provenance.
		
		guardrail_mode: "Implements about a half-dozen safeguards against extremely low-quality samples running for abnormally long times."
		# Previously guardrail_mode set TBProfiler's min % masked to 10% and TBProfiler's min depth to 3, but now these use (100 - sample_max_pct_masked)
		# and sample_min_avg_depth instead.

		output_sample_name: "Override ALL sample names with this string instead."
		# Currently required to deal with certain CDPH edge cases. For Terra data tables, set this to the sample's entity_id column (the one on the far
		# left that acts like an index).
		# TODO: This makes the multi-sample-per-workflow case give the same output. All WDL executers I'm aware of put all scattered outs in different folders,
		# so this won't overwrite per say, but it's not ideal... I need to make sure paired_fastqs being Array[Array[File]] will zip() nicely with an
		# Array[String] version of output_sample_name.

		paired_fastq_sets: "Nested array of paired fastqs, each inner array representing one samples worth of paired fastqs"
		# On a sample-indexed data table on Terra, you'll probably want something like `[[this.read1, this.read2]]`, but you also have the option of running
		# multiple samples at once thanks to the nesting. For example, if you have three samples:
		# [["gs://bucket/foo_1.fq", "gs://bucket/foo_2.fq"], ["gs://bucket/bar_1.fq", "gs://bucket/bar_2.fq"], ["gs://bucket/bizz_1.fq", "gs://bucket/bizz_2.fq"]]
		# Will be read processed as three different samples using WDL scatter(). It is okay if some samples fail and others pass; the passing samples will complete
		# the rest of the pipeline (ie generate diff files) even if other samples "drop out" earlier. This requires some special handling in WDL 1.0 however; please
		# be cautious if you update this pipeline to WDL 1.1 as the workarounds I use for this may require some tinkering under WDL 1.1 standards.

		sample_max_pct_masked: "Samples who have more than this percent (as int, 50 = 50%) of positions with coverage below site_min_depth will be discarded"
		# It'd be more accurate to call this sample_pct_below_site_min_depth but that's far too long!

		sample_min_avg_depth: "If covstats or TBProfiler thinks mean (not median!) depth is below this, throw out this sample"
		# CDC minimum: 50x estimated from read length. We don't estimate from read length, we measure directly, so our values will always be less.

		sample_min_pct_mapped: "If covstats or TBProfiler thinks less than this percent of your sample (after decontam) maps to H37Rv, throw out this sample"
		# This is downstream of fastp cleaning too
		# CDC minimum: 90% mapped to MTBC (not H37Rv) using something like Kraken

		sample_min_q30: "Decontaminated samples with less than this percent (as int, 50 = 50%) of reads above qual score of 30 will be discarded."
		# This runs after decontamination, which by default runs after an initial fastp clean, and is therefore sort of influenced by fastp_avg_qual
		# CDC minimum: 85%. After discussing with CDPH, we're defaulting to 80% instead, because 85% would result in removing so many mostly-good
		# samples that it'd meaningfully affect TB-D's utility for tracking disease.

		site_min_depth: "Positions with coverage below this value will be masked in diff files; see also sample_max_pct_masked"
		# This is explict masking with -, as opposed to "masking to reference"

		skip_covstats: "Should we skip covstats entirely?"
		# Covstats might be entirely removed in a future version as the current version of TBProfiler replaces our old use cases for covstats.

		subsample_cutoff: "If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable, overridden to 500000 (500 GB) if guardrail_mode=True)"
		# One thing SRA taught me is that there will always be somebody who, with the best of intentions, uploads a terabyte of reads as a single sample
	}

	# Flip some QC stuff around
	Int   sample_max_pct_unmapped = 100 - sample_min_pct_mapped
	Float sample_max_pct_masked_float = sample_max_pct_masked / 100.0

	# Some variables we no longer have adjustable by the user to reduce the amount of variable spam on Terra's workflow page
	Boolean QC_soft_pct_mapped = false
	Int quick_tasks_disk_size  = 10  # disk size in GB for file-processing tasks; if running one-workflow-many-samples,
	                                 # increasing this might speed up file localization if you have >5,000 samples
	Boolean strip_all_underscores = false  # currently unused proposed workaround for irregular sample names; 
	                                       # setting this to true can mess with R1/R2 detection if not defined(output_sample_name)
	Boolean tbprofiler_on_bam = just_like_2024
	
	# Used for some workarounds
	String pass = "PASS"

	scatter(paired_fastqs in paired_fastq_sets) {
		call clckwrk_combonation.clean_and_decontam_and_check as decontam_each_sample {
			input:
				CDC_decontamination_reference = false,  # previously decontam_use_CDC_varpipe_ref but not currently supported
				oldschool_docker = just_like_2024,
				unsorted_sam = true,
				force_rename_out = output_sample_name,
				reads_files = paired_fastqs,
				fastp_clean_avg_qual = fastp_avg_qual,
				QC_min_q30 = sample_min_q30,
				strip_all_underscores = strip_all_underscores,
				preliminary_min_q30 = if guardrail_mode then 20 else 1,
				subsample_cutoff = if guardrail_mode then 500000 else subsample_cutoff,
				timeout_map_reads = if guardrail_mode then 300 else 0,
				timeout_decontam = if guardrail_mode then 600 else 0,
				addldisk = if low_resource_mode then 10 else 100,
				memory = if low_resource_mode then 8 else 32
		}

		if(defined(decontam_each_sample.decontaminated_fastq_1)) {
			# This region only executes if decontaminated fastqs exist. We can use this to coerce File? into File by using
			# select_first() where the first element is the File? we know must exist, and the second element never gets used.
			#
			# But what happens if Cromwell changes the behavior of defined() and select_first() so that an undefined File?
			# triggers defined() but not select_first()? This should never happen but given inconsistent handling of optional
			# types in Cromwell over the years, it is worth accounting for.
			#
			# We can't make the second element of the select_first() array completely bogus, because for some reason, all members
			# of select_first() are checked for "is this even remotely valid" even if the first (0 index) member of the
			# select_first() array is defined (eg, is selected). In other words this isn't like a boolean or where you can
			# make the second value some nonsense. (In other words WDL doesn't do short-circuit evaulation for select_first().)
			#
			# So we need something "valid" but still wrong enough to throw an error. The easiest way to do this is to use
			# counts_out_tsv, a File? that is absolutely not going to be valid for tbprofilerFQ. This means Cromwell changing
			# behavior will still bug out, but it will bug out in a way that is immediately detectable.
			#
			File real_decontaminated_fastq_1=select_first([decontam_each_sample.decontaminated_fastq_1, decontam_each_sample.counts_out_tsv]) #!SelectArray
			File real_decontaminated_fastq_2=select_first([decontam_each_sample.decontaminated_fastq_2, decontam_each_sample.counts_out_tsv])
    		
			call tbprofilerFQ_WF.TheiagenTBProfiler as theiagenTBprofilerFQ {
				input:
					fastq1 = real_decontaminated_fastq_1,
					fastq2 = real_decontaminated_fastq_2,
					soft_pct_mapped = false,
					soft_depth = false,
					minimum_median_depth = if guardrail_mode then 3 else 0, # super low to avoid conflict with avg depth
					minimum_mean_depth = sample_min_avg_depth,
					minimum_pct_mapped = sample_min_pct_mapped,
					sample = decontam_each_sample.sample
			}
			# if this sample passes...
			if(theiagenTBprofilerFQ.status_code == pass) {
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
		if(!skip_covstats) {
	
			# covstats to check coverage and percent mapped to reference
			call goleft.covstats as covstats {
				input:
					inputBamOrCram = vcfs_and_bams.left[0],
					allInputIndexes = [vcfs_and_bams.left[1]]
			}
			
			if((covstats.percentUnmapped < sample_max_pct_unmapped) || QC_soft_pct_mapped) {
				if(covstats.coverage > sample_min_avg_depth) {
					
					# make diff files
					call diff.make_mask_and_diff as make_mask_and_diff_after_covstats {
						input:
							bam = vcfs_and_bams.left[0],
							vcf = vcfs_and_bams.right,
							min_coverage_per_site = site_min_depth,
							tbmf = call_as_reference_bedfile,
							max_ratio_low_coverage_sites_per_sample = sample_max_pct_masked_float
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
					min_coverage_per_site = site_min_depth,
					tbmf = call_as_reference_bedfile,
					max_ratio_low_coverage_sites_per_sample = sample_max_pct_masked_float
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

	# NOTE -- currently none of the bam-flavored TBProfiler is a workflow-level out, but I'm leaving this section here
	# in case someone wants to use bam-flavored TBProfiler in the future

	# 1. Coerce bam-flavored TBProfiler into required types
	Array[String] coerced_bam_strains=select_all(profile_bam.sample_and_strain)
	Array[String] coerced_bam_resistances=select_all(profile_bam.sample_and_resistance)
	Array[String] coerced_bam_depths=select_all(profile_bam.sample_and_median_depth)
	
	# 2. Determine if we ran TBProfiler on bams, without relying on defined()
	if(!(length(coerced_bam_strains) == 0)) {
	
		# 3. Determine if we are running on one sample or multiple samples
		if(length(paired_fastq_sets) != 1) {
	
			call sranwrp_processing.cat_strings as collate_bam_strains {      #!UnusedCall
				input:
					strings = coerced_bam_strains,
					out = "strain_reports.txt",
					disk_size = quick_tasks_disk_size
			}
			
			call sranwrp_processing.cat_strings as collate_bam_resistance {   #!UnusedCall
				input:
					strings = coerced_bam_resistances,
					out = "resistance_reports.txt",
					disk_size = quick_tasks_disk_size
			}
	
			call sranwrp_processing.cat_strings as collate_bam_depth {        #!UnusedCall
				input:
					strings = coerced_bam_depths,
					out = "depth_reports.txt",
					disk_size = quick_tasks_disk_size
			}
		}
		
		# if there is only one sample, there's no need to run tasks
		# currently not output and unusued, but I'm leaving them here in case someone needs it later
		if(length(paired_fastq_sets) == 1) {
			String single_sample_tbprof_bam_depth      = coerced_bam_depths[0]      #!UnusedDeclaration
			String single_sample_tbprof_bam_resistance = coerced_bam_resistances[0] #!UnusedDeclaration
			String single_sample_tbprof_bam_strain     = coerced_bam_strains[0]     #!UnusedDeclaration
		}
	}
  	  	
  	# 1. Coerce FQ-flavored TBProfiler into required types
	Array[String] coerced_fq_strains=select_all(theiagenTBprofilerFQ.sample_and_strain)
	Array[String] coerced_fq_resistances=select_all(theiagenTBprofilerFQ.sample_and_resistance)
	Array[String] coerced_fq_depths=select_all(theiagenTBprofilerFQ.sample_and_depth)
	
	# 2. Determine if we ran TBProfiler on FQs, without relying on defined()
	if(!(length(coerced_fq_strains) == 0)) {
	
		# 3. Determine if we are running on one sample or multiple samples
		if(length(paired_fastq_sets) != 1) {

			call sranwrp_processing.cat_strings as collate_fq_strains {     #!UnusedCall
				input:
					strings = coerced_fq_strains,
					out = "strain_reports.txt",
					disk_size = quick_tasks_disk_size
			}
			
			call sranwrp_processing.cat_strings as collate_fq_resistance {  #!UnusedCall
				input:
					strings = coerced_fq_resistances,
					out = "resistance_reports.txt",
					disk_size = quick_tasks_disk_size
			}
			
			call sranwrp_processing.cat_strings as collate_fq_depth {      #!UnusedCall
				input:
					strings = coerced_fq_depths,
					out = "depth_reports.txt",
					disk_size = quick_tasks_disk_size
			}
		}
	
		# if there is only one sample, there's no need to run tasks
		if(length(paired_fastq_sets) == 1) {
			String single_sample_tbprof_fq_depth      = coerced_fq_depths[0]      #!UnusedDeclaration
			String single_sample_tbprof_fq_resistance = coerced_fq_resistances[0] #!UnusedDeclaration
			String single_sample_tbprof_fq_strain     = coerced_fq_strains[0]     #!UnusedDeclaration
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
					
					if((percentUnmapped > sample_max_pct_unmapped) && !(QC_soft_pct_mapped)) { 
						String too_many_unmapped = "COVSTATS_${percentUnmapped}_UNMAPPED_(MAX_${sample_max_pct_unmapped})"
						if(meanCoverage < sample_min_avg_depth) {
							String double_bad = "COVSTATS_BOTH_${percentUnmapped}_UNMAPPED_(MAX_${sample_max_pct_unmapped})_AND_${meanCoverage}_MEAN_COVERAGE_(MIN_${sample_min_avg_depth})"
						} 
					}
					if(meanCoverage < sample_min_avg_depth) {
						String too_low_coverage = "COVSTATS_${meanCoverage}_MEAN_COVERAGE_(MIN_${sample_min_avg_depth})"
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
	# Failed to evaluate 'warnings' (reason 1 of 1): Evaluating flatten(flatten([[select_all(theiagenTBprofilerFQ.warning_codes)], 
	# [select_all(warning_decontam)]])) failed: No coercion defined from wom value(s) '
	# [["EARLYQC_88.112_PCT_ABOVE_Q30_(MIN_0.9)", "EARLYQC_99.61_PCT_MAPPED_(MIN_99.995)"]]' of type 'Array[Array[String]]' to 'Array[String]'.
	#Array[String] warnings = flatten(flatten([[select_all(theiagenTBprofilerFQ.warning_codes)], [select_all(warning_decontam)]]))
		
	output {
		String  tbd_status = select_first([finalcode, pass])
		String? tbd_comment = comment

		# decon/fastp metadata pulled out directly -- only valid if this pipeline ran on a single sample
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
		
		# theiagen!TBProfiler metadata pulled out directly, only valid if single sample
		Float? tbd_qc_median_depth_per_tbprof = theiagenTBprofilerFQ.median_depth[0]
		Float? tbd_qc_avg_depth_per_tbprof = theiagenTBprofilerFQ.avg_depth[0]
		Float? tbd_qc_pct_mapped_per_tbprof = theiagenTBprofilerFQ.pct_reads_mapped[0]
		Float? tbd_qc_pct_genome_covered = theiagenTBprofilerFQ.pct_genome_covered[0]
		Int?    tbd_n_dr_variants = theiagenTBprofilerFQ.n_dr_variants[0]
		Int?    tbd_n_other_variants = theiagenTBprofilerFQ.n_other_variants[0]
		String? tbd_resistance = theiagenTBprofilerFQ.resistance[0]
		String? tbd_strain_per_tbprof = theiagenTBprofilerFQ.strain[0]
		
		# genomic data files
		Array[File]  tbd_bais  = final_bais
		Array[File]  tbd_bams  = final_bams
		Array[File]  tbd_diffs = real_diffs
		Array[File]  tbd_masks = real_masks   # bedgraph
		Array[File]  tbd_vcfs  = minos_vcfs
		
		# metadata files
		Array[File?] tbd_decontam_reports          = decontam_each_sample.counts_out_tsv
		Array[File?] tbd_diff_reports              = real_reports
		Array[File?] tbd_tbprof_fq_jsons           = theiagenTBprofilerFQ.tbprofiler_json
		Array[File?] tbd_tbprof_fq_looker          = theiagenTBprofilerFQ.tbprofiler_looker_csv
		Array[File?] tbd_tbprof_fq_laboratorian    = theiagenTBprofilerFQ.tbprofiler_laboratorian_report_csv
		Array[File?] tbd_tbprof_fq_lims            = theiagenTBprofilerFQ.tbprofiler_lims_report_csv

		# Typically unusued -- these work fine, I just want Terra's UI to be less crowded
		#Array[File?] tbd_tbprof_bam_jsons          = profile_bam.tbprofiler_json
		#Array[File?] tbd_tbprof_bam_summaries      = profile_bam.tbprofiler_txt
		#Array[File?] tbd_covstats_reports          = covstats.covstatsOutfile
		#Float? tbd_qc_pct_unmapped_covstats        = percentUnmapped

		# useful debugging/run information (only valid iff this ran on only one sample)
		String tbd_clockwork_docker      = decontam_each_sample.docker_used[0]
		#String tbd_resistance_coerced   = single_sample_tbprof_fq_resistance
		#String tbd_strain_coerced       = single_sample_tbprof_fq_strain
		#Array[String] pass_or_warnings  = if (length(warnings) > 0) then warnings else ["PASS"] # might not work properly
		#String? tbd_debug_decontam_ERR  = decontam_ERR
		#String? tbd_debug_earlyQC_ERR   = earlyQC_ERR
		#String? tbd_debug_varcall_ERR   = varcall_ERR
		#String? tbd_debug_covstats_ERR  = covstats_ERR
		#String? tbd_debug_vcfdiff_ERR   = vcfdiff_ERR
		#Array[String]? tbd_debug_vcfdiff_errorcode_if_covstats    = vcfdiff_errorcode_if_covstats
		#Array[String]? tbd_debug_vcfdiff_errorcode_if_no_covstats = vcfdiff_errorcode_if_no_covstats
		#Array[String]? tbd_debug_vcfdiff_errorcode_array          = vcfdiff_errorcode_array

	}
}


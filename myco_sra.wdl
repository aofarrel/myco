version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.9.1/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.9.1/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.12/tasks/pull_fastqs.wdl" as sranwrp_pull
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.12/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/tree_nine/0.0.10/tree_nine.wdl" as build_treesWF
import "https://raw.githubusercontent.com/aofarrel/parsevcf/1.1.9/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/fastqc-wdl/0.0.2/fastqc.wdl" as fastqc
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.2.2/tbprofiler_tasks.wdl" as profiler
import "https://raw.githubusercontent.com/aofarrel/TBfastProfiler/main/TBfastProfiler.wdl" as earlyQC


workflow myco {
	input {
		File biosample_accessions
		
		# crash the entire pipeline if any sample is being ornery
		Boolean crash_pipeline_if_varcalling_errors = false
		Boolean crash_pipeline_if_varcalling_timesout = false
		
		# creation + masking of diff files
		Boolean diff_force                   = false
		Int     diff_min_coverage_per_site   = 10
		File?   diff_mask_these_regions
		
		# QC
		Boolean fastqc_on_timeout       = false
		Boolean early_qc_apply_cutoffs  = false
		Float   early_qc_cutoff_q30     = 0.90
		Boolean early_qc_skip_entirely  = false
		
		# shrink large samples
		Int     subsample_cutoff        =  450
		Int     subsample_seed          = 1965
		
		# skip samples that take a very long time
		Int     timeout_decontam_part1  =   20
		Int     timeout_decontam_part2  =   15
		Int     timeout_variant_caller  =  120
		
		# tbprofiler
		Boolean tbprofiler_on_bam       = true
		
		# phylogenetics
		Boolean tree_decoration         = false
		File?   tree_to_decorate
		Float   tree_max_low_coverage_sites = 0.05
		# removed ref genome for tree building
		# eventually remove tree_decoration, will be redundant with tree_to_decorate
		
		# variant caller specifics
		Int  variantcalling_addl_disk    = 100
		Int  variantcalling_cpu          =  16
		Int? variantcalling_mem_height
		Int  variantcalling_memory       =  32
		Int  variantcalling_preemptibles =   1
		Int  variantcalling_retries      =   1
		
	}

	parameter_meta {
		biosample_accessions: "File of BioSample accessions to pull, one accession per line"
		tree_max_low_coverage_sites: "If a diff file has higher than this percent (as float, eg 0.5 = 50%) bad data, do not include it in the tree"
		tree_decoration: "Should usher, taxonium, and NextStrain trees be generated?"
		
		early_qc_apply_cutoffs: "If true, run fastp + TBProfiler on decontaminated fastqs and apply cutoffs to determine which samples should be thrown out."
		early_qc_cutoff_q30: "Decontaminated samples with less than this percentage (as float, 0.5 = 50%) of reads above qual score of 30 will be discarded iff early_qc_apply_cutoffs is also true."
		early_qc_skip_entirely: "Do not run early QC (fastp + fastq-TBProfiler) at all. Does not affect whether or not TBProfiler is later run on bams. Overrides early_qc_apply_cutoffs."
		fastqc_on_timeout: "If true, fastqc one read from a sample when decontamination or variant calling times out"
		
		diff_force: "If true and if tree_decoration is false, generate diff files. (Diff files will always be created if tree_decoration is true.)"
		input_tree: "Base tree to use if tree_decoration = true"
		diff_min_coverage_per_site: "Positions with coverage below this value will be masked in diff files"
		subsample_cutoff: "If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)"
		subsample_seed: "Seed used for subsampling with seqtk"
		tbprofiler_on_bam: "If true, run TBProfiler on BAMs."
		timeout_decontam_part1: "Discard any sample that is still running in clockwork map_reads after this many minutes (set to 0 to never timeout)"
		timeout_decontam_part2: "Discard any sample that is still running in clockwork rm_contam after this many minutes (set to 0 to never timeout)"
		timeout_variant_caller: "Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to 0 to never timeout)"
		diff_mask_these_regions: "Bed file of regions to mask when making diff files"
	}

	# WDL doesn't understand mutual exclusivity, so we have to get a little creative on 
	# our determination of whether or not we want to create diff files.
	if(tree_decoration)  {  Boolean create_diff_files_   = true  }
	if(!tree_decoration) {
		if(!diff_force){  Boolean create_diff_files__  = false }
		if(diff_force) {  Boolean create_diff_files___ = true  }
	}
	Boolean create_diff_files = select_first([create_diff_files_,
											  create_diff_files__, 
											  create_diff_files___])

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
				subsample_cutoff = subsample_cutoff,
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
			out = "pull_reports.txt"
	}

	Array[Array[File]] pulled_fastqs = select_all(paired_fastqs)
	scatter(pulled_fastq in pulled_fastqs) {
		call clckwrk_combonation.combined_decontamination_single_ref_included as decontam_each_sample {
			input:
				unsorted_sam = true,
				reads_files = pulled_fastq,
				timeout_map_reads = timeout_decontam_part1,
				timeout_decontam = timeout_decontam_part2
		}

		if(defined(decontam_each_sample.decontaminated_fastq_1)) {
			# This region only executes if decontaminated fastqs exist. We can use this to coerce File? into File by using
			# select_first() where the first element is the File? we know must exist, and the second element is bogus.
    		File real_decontaminated_fastq_1=select_first([decontam_each_sample.decontaminated_fastq_1, biosample_accessions])
    		File real_decontaminated_fastq_2=select_first([decontam_each_sample.decontaminated_fastq_2, biosample_accessions])

			# if we are NOT skipping earlyQC
			if(!early_qc_skip_entirely) {
				call earlyQC.TBfastProfiler as check_fastqs {
					input:
						fastq1 = real_decontaminated_fastq_1,
						fastq2 = real_decontaminated_fastq_2,
						q30_cutoff = early_qc_cutoff_q30
				}
				
				# if we are filtering out samples via earlyQC...
				if(early_qc_apply_cutoffs) {
					if(check_fastqs.did_this_sample_pass) {
						File possibly_fastp_cleaned_fastq1_passed=select_first([check_fastqs.cleaned_fastq1, real_decontaminated_fastq_1])
				    	File possibly_fastp_cleaned_fastq2_passed=select_first([check_fastqs.cleaned_fastq2, real_decontaminated_fastq_2])
				    	
						call clckwrk_var_call.variant_call_one_sample_ref_included as variant_call_after_earlyQC_filtering {
							input:
								reads_files = [possibly_fastp_cleaned_fastq1_passed, possibly_fastp_cleaned_fastq2_passed],
								
								addldisk = variantcalling_addl_disk,
								cpu = variantcalling_cpu,
								crash_on_error = crash_pipeline_if_varcalling_errors,
								crash_on_timeout = crash_pipeline_if_varcalling_timesout,
								mem_height = variantcalling_mem_height,
								memory = variantcalling_memory,
								preempt = variantcalling_preemptibles,
								retries = variantcalling_retries,
								timeout = timeout_variant_caller
						}
					}
				}
				
				# if we are not filtering out samples via the early qc step (but ran earlyQC anyway)...
				if(!early_qc_apply_cutoffs) {
					File possibly_fastp_cleaned_fastq1=select_first([check_fastqs.cleaned_fastq1, real_decontaminated_fastq_1])
			    	File possibly_fastp_cleaned_fastq2=select_first([check_fastqs.cleaned_fastq2, real_decontaminated_fastq_2])
				    	
					call clckwrk_var_call.variant_call_one_sample_ref_included as variant_call_after_earlyQC_but_not_filtering_samples {
						input:
							reads_files = [possibly_fastp_cleaned_fastq1, possibly_fastp_cleaned_fastq2],
							
							addldisk = variantcalling_addl_disk,
							cpu = variantcalling_cpu,
							crash_on_error = crash_pipeline_if_varcalling_errors,
							crash_on_timeout = crash_pipeline_if_varcalling_timesout,
							mem_height = variantcalling_mem_height,
							memory = variantcalling_memory,
							preempt = variantcalling_preemptibles,
							retries = variantcalling_retries,
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
						crash_on_error = crash_pipeline_if_varcalling_errors,
						crash_on_timeout = crash_pipeline_if_varcalling_timesout,
						mem_height = variantcalling_mem_height,
						memory = variantcalling_memory,
						preempt = variantcalling_preemptibles,
						retries = variantcalling_retries,
						timeout = timeout_variant_caller
				}
			}
		}
			
	}

	# do some wizardry to deal with optionals
	Array[File] minos_vcfs_if_earlyQC_filtered = select_all(variant_call_after_earlyQC_filtering.vcf_final_call_set)
	Array[File] minos_vcfs_if_earlyQC_but_not_filtering = select_all(variant_call_after_earlyQC_but_not_filtering_samples.vcf_final_call_set)
	Array[File] minos_vcfs_if_no_earlyQC = select_all(variant_call_without_earlyQC.vcf_final_call_set)
	
	Array[File] bams_if_earlyQC_filtered = select_all(variant_call_after_earlyQC_filtering.mapped_to_ref)
	Array[File] bams_if_earlyQC_but_not_filtering = select_all(variant_call_after_earlyQC_but_not_filtering_samples.mapped_to_ref)
	Array[File] bams_if_no_earlyQC = select_all(variant_call_without_earlyQC.mapped_to_ref)
	
	Array[File] minos_vcfs = flatten([minos_vcfs_if_earlyQC_filtered, minos_vcfs_if_earlyQC_but_not_filtering, minos_vcfs_if_no_earlyQC])
	Array[File] bams_to_ref = flatten([bams_if_earlyQC_filtered, bams_if_earlyQC_but_not_filtering, bams_if_no_earlyQC])


	scatter(vcfs_and_bams in zip(bams_to_ref, minos_vcfs)) {
		call diff.make_mask_and_diff as make_mask_and_diff {
			input:
				bam = vcfs_and_bams.left,
				vcf = vcfs_and_bams.right,
				min_coverage_per_site = diff_min_coverage_per_site,
				tbmf = diff_mask_these_regions,
				diffs = create_diff_files
		}
		
		if(tbprofiler_on_bam) {
			call profiler.tb_profiler_bam as profile_bam {
					input:
						bam = vcfs_and_bams.left
			}
		}
	}

	if(defined(profile_bam.strain)) {
		Array[String] coerced_strains=select_all(profile_bam.strain)
		Array[String] coerced_resistance=select_all(profile_bam.resistance)
		Array[String] coerced_depth=select_all(profile_bam.median_depth)

		call sranwrp_processing.cat_strings as collate_bam_strains {
			input:
				strings = coerced_strains,
				out = "strain_reports.txt"
		}
		
		call sranwrp_processing.cat_strings as collate_bam_resistance {
			input:
				strings = coerced_resistance,
				out = "resistance_reports.txt"
		}

		call sranwrp_processing.cat_strings as collate_bam_depth {
			input:
				strings = coerced_depth,
				out = "depth_reports.txt"
		}
  	}
	
	

	if(fastqc_on_timeout) {
		Array[File] bad_fastqs_decontam_ = select_all(decontam_each_sample.check_this_fastq)
		Array[File] bad_fastqs_varcallr_ = select_all(flatten([variant_call_without_earlyQC.check_this_fastq, variant_call_after_earlyQC_but_not_filtering_samples.check_this_fastq]))
		Array[Array[File]] bad_fastqs_   = [bad_fastqs_decontam_, bad_fastqs_varcallr_]
		if(length(decontam_each_sample.check_this_fastq)>=1 && length(bad_fastqs_varcallr_)>=1) {
			Array[File] bad_fastqs_both  = flatten(bad_fastqs_)  
		}
		if(length(decontam_each_sample.check_this_fastq)>=1) {
			Array[File] bad_fastqs_decontam = select_all(bad_fastqs_decontam_)
		}
		if(length(bad_fastqs_varcallr_)>=1) {
			Array[File] bad_fastqs_varcallr = select_all(bad_fastqs_varcallr_)
		}
		Array[File] fastqs = select_first([bad_fastqs_both, bad_fastqs_decontam, bad_fastqs_varcallr])
		if(length(fastqs)>0) {
			call fastqc.FastqcWF {
				input:
					fastqs = fastqs
			}
		}
		
	}

	if(tree_decoration) {
		# diff files must exist if tree_decoration is true, so we can force the Array[File?]?
		# into an Array[File] with the classic "select_first() with a bogus fallback" hack
		Array[File] coerced_diffs = select_first([select_all(make_mask_and_diff.diff), minos_vcfs])
		Array[File] coerced_reports = select_first([select_all(make_mask_and_diff.report), minos_vcfs])
		call build_treesWF.Tree_Nine as trees {
			input:
				diffs = coerced_diffs,
				input_tree = tree_to_decorate,
				coverage_reports = coerced_reports,
				max_low_coverage_sites = tree_max_low_coverage_sites
		}
	}

	output {
		# raw files
		Array[File?] diffs = make_mask_and_diff.diff
		Array[File]  masks = make_mask_and_diff.mask_file
		Array[File]  vcfs  = minos_vcfs
		
		# metadata
		File          download_report        = merge_reports.outfile
		Array[File]?  fastqc_reports         = FastqcWF.reports
		Array[File?]? fastp_reports          = check_fastqs.fastp_txt
		File?         tbprof_bam_depths      = collate_bam_depth.outfile
		Array[File?]? tbprof_bam_jsons       = profile_bam.tbprofiler_json
		File?         tbprof_bam_strains     = collate_bam_strains.outfile
		Array[File?]? tbprof_bam_summaries   = profile_bam.tbprofiler_txt
		File?         tbprof_bam_resistances = collate_bam_resistance.outfile
		
		# tree nine
		File?        tree_nwk         = trees.tree_nwk
		File?        tree_usher       = trees.tree_usher_raw
		File?        tree_taxonium    = trees.tree_taxonium
		File?        tree_nextstrain  = trees.tree_nextstrain
		Array[File]? trees_nextstrain = trees.subtrees_nextstrain
	}
}
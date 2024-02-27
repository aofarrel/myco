version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.12.2/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.12.2/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.18/tasks/pull_fastqs.wdl" as sranwrp_pull
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.18/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/tree_nine/0.0.15/tree_nine.wdl" as build_treesWF
import "https://raw.githubusercontent.com/aofarrel/vcf_to_diff_wdl/0.0.3/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.2.3/tbprofiler_tasks.wdl" as profiler
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.2.4/thiagen_tbprofiler.wdl" as tbprofilerFQ_WF # fka earlyQC

import "https://raw.githubusercontent.com/aofarrel/TBfastProfiler/0.0.6/TBfastProfiler.wdl" as qc_fastqsWF # aka earlyQC
import "https://raw.githubusercontent.com/aofarrel/goleft-wdl/0.1.2/goleft_functions.wdl" as goleft


workflow myco {
	input {
		File biosample_accessions
	
		# covstats (occurs after variant calling)
		Float   covstatsQC_max_percent_unmapped = 10
		Float   covstatsQC_minimum_coverage = 2
		Boolean covstatsQC_skip_entirely = false
		
		# creation + masking of diff files
		Int     diffQC_low_coverage_cutoff      = 10
		File?   diffQC_mask_bedfile
		Int     diffQC_max_percent_low_coverage =  5
		
		# QC
		Boolean earlyQC_skip_QC             = false
		Int     earlyQC_minimum_percent_q30 = 90
		Boolean earlyQC_skip_trimming       = false
		Int     earlyQC_trim_qual_below     = 30
		Boolean earlyQC_skip_entirely       = true
		
		# shrink large samples
		Int     subsample_cutoff        =  450
		Int     subsample_seed          = 1965
		
		# skip samples that take a very long time
		Int     timeout_decontam_part1  =   20
		Int     timeout_decontam_part2  =   15
		Int     timeout_variant_caller  =  120
		
		# tbprofiler on BAMs
		Boolean tbprofiler_on_bam       = true
		
		# phylogenetics
		Boolean tree_decoration         = false
		File?   tree_to_decorate

		# runtime attributes
		Int quick_tasks_disk_size               =  10
		Int  variantcalling_addl_disk           = 100
		Int  variantcalling_cpu                 =  16
		Boolean variantcalling_crash_on_error   = false
		Boolean variantcalling_crash_on_timeout = false
		Int? variantcalling_mem_height
		Int  variantcalling_memory              =  32
		Int  variantcalling_preemptibles        =   1
		Int  variantcalling_retries             =   1
		Boolean variantcalling_ssd              = true
		
	}

	parameter_meta {
		biosample_accessions: "File of BioSample accessions to pull, one accession per line"
		covstatsQC_minimum_coverage: "If covstats thinks MEAN coverage is below this, throw out this sample - not to be confused with TBProfiler MEDIAN coverage"
		covstatsQC_max_percent_unmapped: "If covstats thinks this percent (as int, 50 = 50%) of data does not map to H37Rv, throw out this sample"
		covstatsQC_skip_entirely: "Should we skip covstats entirely?"
		diffQC_mask_bedfile: "Bed file of regions to mask when making diff files (default: R00000039_repregions.bed)"
		diffQC_max_percent_low_coverage: "Samples who have more than this percent (as int, 50 = 50%) of positions with coverage below diffQC_low_coverage_cutoff will be discarded"
		diffQC_low_coverage_cutoff: "Positions with coverage below this value will be masked in diff files"
		earlyQC_minimum_percent_q30: "Decontaminated samples with less than this percent (as int, 50 = 50%) of reads above qual score of 30 will be discarded. Negated by earlyQC_skip_QC or earlyQC_skip_entirely being false."
		earlyQC_skip_entirely: "Do not run earlyQC (fastp + fastq-TBProfiler) at all - no trimming, no QC, nothing. Does not affect tbprofiler_on_bam."
		earlyQC_skip_QC: "Run earlyQC (unless earlyQC_skip_entirely is true), but do not throw out samples that fail QC. Independent of earlyQC_skip_trimming."
		earlyQC_skip_trimming: "Run earlyQC (unless earlyQC_skip_entirely is true), and remove samples that fail QC (unless earlyQC_skip_QC is true), but do not use fastp's cleaned fastqs."
		earlyQC_trim_qual_below: "Trim reads with an average quality score below this value. Independent of earlyQC_minimum_percent_q30. Overridden by earlyQC_skip_trimming or earlyQC_skip_entirely being true."
		quick_tasks_disk_size: "Disk size in GB to use for quick file-processing tasks; increasing this might slightly speed up file localization"
		subsample_cutoff: "If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)"
		subsample_seed: "Seed used for subsampling with seqtk"
		tbprofiler_on_bam: "If true, run TBProfiler on BAMs"
		timeout_decontam_part1: "Discard any sample that is still running in clockwork map_reads after this many minutes (set to 0 to never timeout)"
		timeout_decontam_part2: "Discard any sample that is still running in clockwork rm_contam after this many minutes (set to 0 to never timeout)"
		timeout_variant_caller: "Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to 0 to never timeout)"
		tree_decoration: "Should usher, taxonium, and NextStrain trees be generated?"
		tree_to_decorate: "Base tree to use if tree_decoration = true"		
	}
	
	# convert percent integers to floats (except covstatsQC_max_percent_unmapped)
	Float diffQC_max_percent_low_coverage_float = diffQC_max_percent_low_coverage / 100
	Float earlyQC_minimum_percent_q30_float = earlyQC_minimum_percent_q30 / 100

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
			out = "pull_reports.txt",
			disk_size = quick_tasks_disk_size
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
			if(!earlyQC_skip_entirely) {
				call qc_fastqsWF.TBfastProfiler as qc_fastqs {
					input:
						fastq1 = real_decontaminated_fastq_1,
						fastq2 = real_decontaminated_fastq_2,
						q30_cutoff = earlyQC_minimum_percent_q30_float,
						average_qual = earlyQC_trim_qual_below,
						output_fastps_cleaned_fastqs = !(earlyQC_skip_trimming)
				}
				
				# if we are filtering out samples via earlyQC...
				if(!(earlyQC_skip_QC)) {
					if(qc_fastqs.did_this_sample_pass) {
						File possibly_fastp_cleaned_fastq1_passed=select_first([qc_fastqs.cleaned_fastq1, real_decontaminated_fastq_1])
				    	File possibly_fastp_cleaned_fastq2_passed=select_first([qc_fastqs.cleaned_fastq2, real_decontaminated_fastq_2])
				    	
						call clckwrk_var_call.variant_call_one_sample_ref_included as variant_call_after_earlyQC_filtering {
							input:
								reads_files = [possibly_fastp_cleaned_fastq1_passed, possibly_fastp_cleaned_fastq2_passed],
								
								addldisk = variantcalling_addl_disk,
								cpu = variantcalling_cpu,
								crash_on_error = variantcalling_crash_on_error,
								crash_on_timeout = variantcalling_crash_on_timeout,
								mem_height = variantcalling_mem_height,
								memory = variantcalling_memory,
								preempt = variantcalling_preemptibles,
								retries = variantcalling_retries,
								ssd = variantcalling_ssd,
								timeout = timeout_variant_caller
						}
					}
				}
				
				# if we are not filtering out samples via the early qc step (but ran earlyQC anyway)...
				if(earlyQC_skip_QC) {
					File possibly_fastp_cleaned_fastq1=select_first([qc_fastqs.cleaned_fastq1, real_decontaminated_fastq_1])
			    	File possibly_fastp_cleaned_fastq2=select_first([qc_fastqs.cleaned_fastq2, real_decontaminated_fastq_2])
				    	
					call clckwrk_var_call.variant_call_one_sample_ref_included as variant_call_after_earlyQC_but_not_filtering_samples {
						input:
							reads_files = [possibly_fastp_cleaned_fastq1, possibly_fastp_cleaned_fastq2],
							
							addldisk = variantcalling_addl_disk,
							cpu = variantcalling_cpu,
							crash_on_error = variantcalling_crash_on_error,
							crash_on_timeout = variantcalling_crash_on_timeout,
							mem_height = variantcalling_mem_height,
							memory = variantcalling_memory,
							preempt = variantcalling_preemptibles,
							retries = variantcalling_retries,
							ssd = variantcalling_ssd,
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
						mem_height = variantcalling_mem_height,
						memory = variantcalling_memory,
						preempt = variantcalling_preemptibles,
						retries = variantcalling_retries,
						ssd = variantcalling_ssd,
						timeout = timeout_variant_caller
				}
			}
		}
			
	}

	# do some wizardry to deal with optionals
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
	
	# Now we need to essentially scatter on three arrays: minos_vcfs, bams, and bais. This is trivial in CWL, but
	# isn't intutive in WDL -- in fact, it's arguably not possible!
	#
	# In WDL, you *should* be able to define a custom struct (think of it like a Python object), which can be seen
	# in the WDL spec (https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#struct-definition), and
	# scatter on that struct. However, in my experience, Cromwell gets pretty buggy when you try to scatter on structs
	# especially if those structs contain files (and our structs most definitely contain files). We can't use zip() 
	# either, because that only works if you have two arrays, not three.
	#
	# So, naturally, we're going to do something cheesy.
	
	Array[Array[File]] bams_and_bais = [final_bams, final_bais]
	Array[Array[File]] bam_per_bai = transpose(bams_and_bais)
	
	# bams_and_bais might look like this: [[SAM1234.bam, SAM1235.bam], [SAM1234.bam.bai, SAM1235.bam.bai]]
	# bams_per_bais might look like this: [[SAM1234.bam, SAM1234.bam.bai], [SAM1235.bam, SAM1235.bam.bai]]

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
	
	# TODO: compare select_first() used here to nested flatten(select_all(), select_all()) used in myco_raw -- are there nulls in output?
	Array[File?] real_diffs = select_first([make_mask_and_diff_after_covstats.diff, make_mask_and_diff_no_covstats.diff])
	Array[File?] real_reports = select_first([make_mask_and_diff_after_covstats.report, make_mask_and_diff_no_covstats.report])
	Array[File?] real_masks = select_first([make_mask_and_diff_after_covstats.mask_file, make_mask_and_diff_no_covstats.mask_file])
	
	# these will not crash even if profile_bam/qc_fastqs did not run (see SOTHWO explaination in myco_raw.wdl)
	Array[String] coerced_bam_strains=select_all(profile_bam.strain)
	Array[String] coerced_bam_resistance=select_all(profile_bam.resistance)
	Array[String] coerced_bam_depth=select_all(profile_bam.median_depth)
	Array[String] coerced_fq_strains=select_all(qc_fastqs.samp_strain)
	Array[String] coerced_fq_resistance=select_all(qc_fastqs.samp_resistance)

	# pull TBProfiler information, if we ran TBProfiler on bams
	if(!(length(coerced_bam_strains) == 0)) {
	
		call sranwrp_processing.cat_strings as collate_bam_strains {
			input:
				strings = coerced_bam_strains,
				out = "strain_reports.txt",
				disk_size = quick_tasks_disk_size
		}
		
		call sranwrp_processing.cat_strings as collate_bam_resistance {
			input:
				strings = coerced_bam_resistance,
				out = "resistance_reports.txt",
				disk_size = quick_tasks_disk_size
		}

		call sranwrp_processing.cat_strings as collate_bam_depth {
			input:
				strings = coerced_bam_depth,
				out = "depth_reports.txt",
				disk_size = quick_tasks_disk_size
		}
  	}
  	
  	# pull TBProfiler information, if we ran TBProfiler on fastqs
  	if(!(length(coerced_fq_strains) == 0)) {

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
		if(length(real_diffs)>0) {
			# diff files must exist if tree_decoration is true (unless no samples passed) 
			# so we can force the Array[File?]? into an Array[File] with the classic 
			# "select_first() with a bogus fallback" hack
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
		File          download_report        = merge_reports.outfile
		Array[File?]  fastp_reports          = qc_fastqs.fastp_txt
		File?         tbprof_bam_depths      = collate_bam_depth.outfile
		Array[File?]  tbprof_bam_jsons       = profile_bam.tbprofiler_json
		File?         tbprof_bam_strains     = collate_bam_strains.outfile
		Array[File?]  tbprof_bam_summaries   = profile_bam.tbprofiler_txt
		File?         tbprof_bam_resistances = collate_bam_resistance.outfile
		Array[File?]  tbprof_fq_jsons        = qc_fastqs.tbprofiler_json
		Array[File?]  tbprof_fq_looker       = qc_fastqs.tbprofiler_looker_csv
		Array[File?]  tbprof_fq_laboratorian = qc_fastqs.tbprofiler_laboratorian_report_csv
		File?         tbprof_fq_strains      = collate_fq_strains.outfile
		Array[File?]  tbprof_fq_summaries    = qc_fastqs.tbprofiler_txt
		File?         tbprof_fq_resistances  = collate_fq_resistance.outfile
		
		# tree nine
		File?        tree_nwk         = trees.tree_nwk
		File?        tree_usher       = trees.tree_usher
		File?        tree_taxonium    = trees.tree_taxonium
		File?        tree_nextstrain  = trees.tree_nextstrain
		Array[File]? trees_nextstrain = trees.subtrees_nextstrain
	}
}

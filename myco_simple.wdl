version 1.0
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.12.1/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/fastp-wdl/main/fastp_tasks.wdl" as fastp
import "https://raw.githubusercontent.com/aofarrel/tree_nine/0.0.12/tree_nine.wdl" as build_treesWF
import "https://raw.githubusercontent.com/aofarrel/vcf_to_diff_wdl/0.0.1/vcf_to_diff.wdl" as diff

# This is a stripped-down version of myco which only runs the variant caller, vcf-to-diff, and (optionally) fastp and/or Tree Nine.
# There is NO decontamination, TB-Profiler, covstats, or QC beyond "does the variant caller crash/time out or not."
# Need to decontaminate your fastqs, but don't want to run myco_raw? Try Decontam_And_Combine_One_Samples_Fastqs, which can be found
# on Dockstore: https://dockstore.org/workflows/github.com/aofarrel/clockwork-wdl/Decontam_And_Combine_One_Samples_Fastqs


workflow myco {
	input {
		Array[Array[File]] paired_decontaminated_fastq_sets

		Boolean decorate_tree      = false
		Boolean fastp_clean        = true
		File?   input_tree
		Float   max_low_coverage_sites = 0.05
		Int     min_coverage_per_site = 10
		Int     timeout_variant_caller =  120
		File?   typical_tb_masked_regions
	}

	parameter_meta {
		decorate_tree: "Should usher, taxonium, and NextStrain trees be generated? Requires input_tree and ref_genome"
		fastp_clean: "If true, clean reads with fastp before calling variants"
		input_tree: "Base tree to use if decorate_tree = true"
		max_low_coverage_sites: "If a diff file has higher than this percent (0.5 = 50%) bad data, do not include it in the tree"
		min_coverage_per_site: "Positions with coverage below this value will be masked in diff files"
		paired_decontaminated_fastq_sets: "Nested array of decontaminated and merged fastq pairs. Each inner array represents one sample; each sample needs precisely one gzipped forward read and one gzipped reverse read."
		subsample_cutoff: "If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)"
		subsample_seed: "Seed used for subsampling with seqtk"
		timeout_variant_caller: "Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to 0 to never timeout)"
		typical_tb_masked_regions: "Bed file of regions to mask when making diff files"
	}

	scatter(paired_fastqs in paired_decontaminated_fastq_sets) {
		if(fastp_clean) {
			call fastp.merge_then_fastp as clean {
				input:
					reads_files = paired_fastqs
			}
		}
		Array[File] cleaned_or_original_fqs = if fastp_clean then select_all([clean.very_clean_fastq1, clean.very_clean_fastq2]) else paired_fastqs
	
		call clckwrk_var_call.variant_call_one_sample_ref_included as variant_call_each_sample {
			input:
				reads_files = cleaned_or_original_fqs,
				timeout = timeout_variant_caller
		}
	}

	Array[File] minos_vcfs=select_all(variant_call_each_sample.adjudicated_vcf)
	Array[File] bams_to_ref=select_all(variant_call_each_sample.bam)


	scatter(vcfs_and_bams in zip(bams_to_ref, minos_vcfs)) {
		call diff.make_mask_and_diff as make_mask_and_diff {
			input:
				bam = vcfs_and_bams.left,
				vcf = vcfs_and_bams.right,
				min_coverage_per_site = min_coverage_per_site,
				max_ratio_low_coverage_sites_per_sample = max_low_coverage_sites,
				tbmf = typical_tb_masked_regions
			}
	}

	if(decorate_tree) {
		if(length(make_mask_and_diff.diff) > 0) {
			Array[File] coerced_diffs = select_all(make_mask_and_diff.diff)
			call build_treesWF.Tree_Nine as trees {
				input:
					diffs = coerced_diffs,
					input_tree = input_tree,
			}
		}
	}

	output {
		Array[File] minos = minos_vcfs
		Array[File] masks = make_mask_and_diff.mask_file
		Array[File?] diffs = make_mask_and_diff.diff

		# tree nine
		File? tree_nwk = trees.tree_nwk
		File? tree_usher = trees.tree_usher_raw
		File? tree_taxonium = trees.tree_taxonium
		File? tree_nextstrain = trees.tree_nextstrain
		Array[File]? trees_nextstrain = trees.subtrees_nextstrain
	}
}

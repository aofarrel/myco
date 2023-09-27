version 1.0
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.9.2/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/tree_nine/0.0.10/tree_nine.wdl" as build_treesWF
import "https://raw.githubusercontent.com/aofarrel/parsevcf/1.1.8/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/fastqc-wdl/0.0.2/fastqc.wdl" as fastqc

# This is a stripped-down version of myco which only runs the variant caller, vcf-to-diff, and (optionally) Tree Nine and/or fastQC.
# There is NO decontamination, fastp, TB-Profiler, covstats, or QC beyond "does the variant caller crash/time out or not."
# Need to decontaminate your fastqs, but don't want to run myco_raw? Try Decontam_And_Combine_One_Samples_Fastqs, which can be found
# on Dockstore: https://dockstore.org/workflows/github.com/aofarrel/clockwork-wdl/Decontam_And_Combine_One_Samples_Fastqs


workflow myco {
	input {
		Array[Array[File]] paired_decontaminated_fastq_sets

		Boolean decorate_tree      = false
		Boolean fastqc_on_timeout  = false
		Boolean force_diff         = false
		File?   input_tree
		Float   max_low_coverage_sites = 0.05
		Int     min_coverage_per_site = 10
		File?   ref_genome_for_tree_building
		Int     timeout_variant_caller =  120
		File?   typical_tb_masked_regions
	}

	parameter_meta {
		decorate_tree: "Should usher, taxonium, and NextStrain trees be generated? Requires input_tree and ref_genome"
		fastqc_on_timeout: "If true, fastqc one read from a sample when decontamination or variant calling times out"
		force_diff: "If true and if decorate_tree is false, generate diff files. (Diff files will always be created if decorate_tree is true.)"
		input_tree: "Base tree to use if decorate_tree = true"
		max_low_coverage_sites: "If a diff file has higher than this percent (0.5 = 50%) bad data, do not include it in the tree"
		min_coverage_per_site: "Positions with coverage below this value will be masked in diff files"
		paired_decontaminated_fastq_sets: "Nested array of decontaminated and merged fastq pairs. Each inner array represents one sample; each sample needs precisely one forward read and one reverse read."
		ref_genome_for_tree_building: "Ref genome for building trees -- must have ONLY `>NC_000962.3` on its first line"
		subsample_cutoff: "If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)"
		subsample_seed: "Seed used for subsampling with seqtk"
		timeout_variant_caller: "Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to 0 to never timeout)"
		typical_tb_masked_regions: "Bed file of regions to mask when making diff files"
	}

	# WDL doesn't understand mutual exclusivity, so we have to get a little creative on 
	# our determination of whether or not we want to create diff files.
	if(decorate_tree)  {  Boolean create_diff_files_   = true  }
	if(!decorate_tree) {
		if(!force_diff){  Boolean create_diff_files__  = false }
		if(force_diff) {  Boolean create_diff_files___ = true  }
	}
	Boolean create_diff_files = select_first([create_diff_files_,
											  create_diff_files__, 
											  create_diff_files___])
    

	scatter(paired_fastqs in paired_decontaminated_fastq_sets) {
			call clckwrk_var_call.variant_call_one_sample_ref_included as variant_call_each_sample {
				input:
					reads_files = paired_fastqs,
					timeout = timeout_variant_caller
			}
		}

	if(fastqc_on_timeout) {
		# Note: This might be problematic in some situations -- may need to make this look like myco_sra
		# But until then, I'm going to stick with this simpler implementation
		if(length(variant_call_each_sample.check_this_fastq)>1) {
			Array[File] bad_fastqs_varcallr = select_all(variant_call_each_sample.check_this_fastq)
		}
		call fastqc.FastqcWF {
			input:
				fastqs = select_first([bad_fastqs_varcallr])
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
				tbmf = typical_tb_masked_regions,
				diffs = create_diff_files
		}
	}

	if(decorate_tree) {
		# diff files must exist if decorate_tree is true, so we can force the Array[File?]?
		# into an Array[File] with the classic "select_first() with a bogus fallback" hack
		Array[File] coerced_diffs = select_first([select_all(make_mask_and_diff.diff), minos_vcfs])
		Array[File] coerced_reports = select_first([select_all(make_mask_and_diff.report), minos_vcfs])
		call build_treesWF.Tree_Nine as trees {
			input:
				diffs = coerced_diffs,
				input_tree = input_tree,
				ref_genome = ref_genome_for_tree_building,
				coverage_reports = coerced_reports,
				max_low_coverage_sites = max_low_coverage_sites
		}
	}

	output {
		Array[File] minos = minos_vcfs
		Array[File] masks = make_mask_and_diff.mask_file
		Array[File?] diffs = make_mask_and_diff.diff
		Array[File]? fastqc_reports = FastqcWF.reports

		# tree nine
		File? tree_nwk = trees.tree_nwk
		File? tree_usher = trees.tree_usher_raw
		File? tree_taxonium = trees.tree_taxonium
		File? tree_nextstrain = trees.tree_nextstrain
		Array[File]? trees_nextstrain = trees.subtrees_nextstrain
	}
}

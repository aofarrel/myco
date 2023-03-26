version 1.0
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.7.0/workflows/refprep-TB.wdl" as clockwork_ref_prepWF
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.7.0/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/tree-nine/0.0.4/usher_sampled.wdl" as build_treesWF
import "https://raw.githubusercontent.com/aofarrel/parsevcf/1.1.4/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/fastqc-wdl/main/fastqc.wdl" as fastqc

workflow myco {
	input {
		Array[Array[File]] paired_decontaminated_fastq_sets
		File typical_tb_masked_regions

		Float   bad_data_threshold = 0.05
		Boolean decorate_tree      = false
		Boolean fastqc_on_timeout  = false
		Boolean force_diff         = false
		File?   input_tree
		Int     min_coverage = 10
		File?   ref_genome_for_tree_building
		Int     subsample_cutoff       =  450
		Int     subsample_seed         = 1965
		Int     timeout_variant_caller =  120
	}

	parameter_meta {
		bad_data_threshold: "If a diff file has higher than this percent (0.5 = 50%) bad data, do not include it in the tree"
		decorate_tree: "Should usher, taxonium, and NextStrain trees be generated? Requires input_tree and ref_genome"
		fastqc_on_timeout: "If true, fastqc one read from a sample when decontamination or variant calling times out"
		force_diff: "If true and if decorate_tree is false, generate diff files. (Diff files will always be created if decorate_tree is true.)"
		input_tree: "Base tree to use if decorate_tree = true"
		min_coverage: "Positions with coverage below this value will be masked in diff files"
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
    
    call clockwork_ref_prepWF.ClockworkRefPrepTB

	scatter(paired_fastqs in paired_decontaminated_fastq_sets) {
			call clckwrk_var_call.variant_call_one_sample_simple as varcall_with_array {
				input:
					ref_dir = ClockworkRefPrepTB.tar_indexd_H37Rv_ref,
					reads_files = paired_fastqs,
					timeout = timeout_variant_caller
			}
		}

	if(fastqc_on_timeout) {
		# Note: This might be problematic in some situations -- may need to make this look like myco_sra
		# But until then, I'm going to stick with this simpler implementation
		if(length(varcall_with_array.check_this_fastq)>1) {
			Array[File] bad_fastqs_varcallr = select_all(varcall_with_array.check_this_fastq)
		}
		call fastqc.FastqcWF {
			input:
				fastqs = select_first([bad_fastqs_varcallr])
		}
	}

	Array[File] minos_vcfs=select_all(varcall_with_array.vcf_final_call_set)
	Array[File] bams_to_ref=select_all(varcall_with_array.mapped_to_ref)


	scatter(vcfs_and_bams in zip(bams_to_ref, minos_vcfs)) {
		call diff.make_mask_and_diff as make_mask_and_diff {
			input:
				bam = vcfs_and_bams.left,
				vcf = vcfs_and_bams.right,
				min_coverage = min_coverage,
				tbmf = typical_tb_masked_regions,
				diffs = create_diff_files
		}
	}

	if(decorate_tree) {
		# diff files must exist if decorate_tree is true, so we can force the Array[File?]?
		# into an Array[File] with the classic "select_first() with a bogus fallback" hack
		Array[File] coerced_diffs = select_first([select_all(make_mask_and_diff.diff), minos_vcfs])
		Array[File] coerced_reports = select_first([select_all(make_mask_and_diff.report), minos_vcfs])
		call build_treesWF.usher_sampled_diff_to_taxonium as trees {
			input:
				diffs = coerced_diffs,
				input_mutation_annotated_tree = input_tree,
				ref = ref_genome_for_tree_building,
				coverage_reports = coerced_reports,
				bad_data_threshold = bad_data_threshold
		}
	}

	output {
		Array[File] minos = minos_vcfs
		Array[File] masks = make_mask_and_diff.mask_file
		Array[File?] diffs = make_mask_and_diff.diff
		File? tree_usher = trees.usher_tree
		File? tree_taxonium = trees.taxonium_tree
		File? tree_nextstrain = trees.nextstrain_tree
		Array[File]? trees_nextstrain = trees.nextstrain_trees
		Array[File]? fastqc_reports = FastqcWF.reports
	}
}

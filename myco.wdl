version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.7.0/workflows/refprep-TB.wdl" as clockwork_ref_prepWF
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.7.0/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.7.0/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/usher-sampled-wdl/0.0.2/usher_sampled.wdl" as build_treesWF
import "https://raw.githubusercontent.com/aofarrel/parsevcf/1.1.2/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/fastqc-wdl/main/fastqc.wdl" as fastqc

workflow myco {
	input {
		Array[Array[File]] paired_fastq_sets
		File typical_tb_masked_regions

		Float   bad_data_threshold = 0.05
		Boolean decorate_tree      = false
		Boolean fastqc_on_timeout  = false
		File?   input_tree
		Int     min_coverage = 10
		File?   ref_genome_for_tree_building
		Int     subsample_cutoff       =  450
		Int     subsample_seed         = 1965
		Int     timeout_decontam_part1 =   20
		Int     timeout_decontam_part2 =   15
		Int     timeout_variant_caller =  120
	}

	parameter_meta {
		bad_data_threshold: "If a diff file has higher than this percent (0.5 = 50%) bad data, do not include it in the tree"
		decorate_tree: "Should usher, taxonium, and NextStrain trees be generated? Requires input_tree and ref_genome"
		fastqc_on_timeout: "If true, fastqc one read from a sample when decontamination times out (see timeout_decontam)"
		input_tree: "Base tree to use if decorate_tree = true"
		min_coverage: "Positions with coverage below this value will be masked in diff files"
		paired_fastq_sets: "Nested array of paired fastqs, each inner array representing one samples worth of paired fastqs"
		ref_genome_for_tree_building: "Ref genome for building trees -- must have ONLY `>NC_000962.3` on its first line"
		subsample_cutoff: "If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)"
		subsample_seed: "Seed used for subsampling with seqtk"
		timeout_decontam_part1: "Discard any sample that is still running in clockwork map_reads after this many minutes (set to -1 to never timeout)"
		timeout_decontam_part2: "Discard any sample that is still running in clockwork rm_contam after this many minutes (set to -1 to never timeout)"
		timeout_variant_caller: "Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to -1 to never timeout)"
		typical_tb_masked_regions: "Bed file of regions to mask when making diff files"
	}

	call clockwork_ref_prepWF.ClockworkRefPrepTB

	scatter(paired_fastqs in paired_fastq_sets) {
		call clckwrk_combonation.combined_decontamination_single as per_sample_decontam {
			input:
				unsorted_sam = true,
				reads_files = paired_fastqs,
				tarball_ref_fasta_and_index = ClockworkRefPrepTB.tar_indexd_dcontm_ref,
				ref_fasta_filename = "ref.fa",
				filename_metadata_tsv = "remove_contam_metadata.tsv",
				subsample_cutoff = subsample_cutoff,
				subsample_seed = subsample_seed,
				timeout_map_reads = timeout_decontam_part1,
				timeout_decontam = timeout_decontam_part2
		}

		if(defined(per_sample_decontam.decontaminated_fastq_1)) {
			# This region only executes if decontaminated fastqs exist.
			# We can use this to coerce File? into File by using a
			# select_first() where the first element is the File? we know
			# absolutely must exist, and the second element is bogus.
			File real_decontaminated_fastq_1=select_first([per_sample_decontam.decontaminated_fastq_1, 
					typical_tb_masked_regions])
			File real_decontaminated_fastq_2=select_first([per_sample_decontam.decontaminated_fastq_2, 
					typical_tb_masked_regions])

			call clckwrk_var_call.variant_call_one_sample_simple as varcall_with_array {
				input:
					ref_dir = ClockworkRefPrepTB.tar_indexd_H37Rv_ref,
					reads_files = [real_decontaminated_fastq_1, real_decontaminated_fastq_2],
					timeout = timeout_variant_caller
			}
		}

	}

	if(fastqc_on_timeout) {
		if(length(per_sample_decontam.check_this_fastq)>1 && length(varcall_with_array.check_this_fastq)>1) {
			Array[File] bad_fastqs_both = select_all(per_sample_decontam.check_this_fastq)
		}
		if(length(per_sample_decontam.check_this_fastq)>1) {
			Array[File] bad_fastqs_decontam = select_all(per_sample_decontam.check_this_fastq)
		}
		if(length(varcall_with_array.check_this_fastq)>1) {
			Array[File] bad_fastqs_varcallr = select_all(varcall_with_array.check_this_fastq)
		}
		call fastqc.FastqcWF {
			input:
				fastqs = select_first([bad_fastqs_both, bad_fastqs_decontam, bad_fastqs_varcallr])
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
				tbmf = typical_tb_masked_regions
		}
	}

	if(decorate_tree) {
		call build_treesWF.usher_sampled_diff_to_taxonium as trees {
			input:
				diffs = make_mask_and_diff.diff,
				i = input_tree,
				ref = ref_genome_for_tree_building,
				coverage_reports = make_mask_and_diff.report,
				bad_data_threshold = bad_data_threshold
		}
	}

	output {
		Array[File] minos = minos_vcfs
		Array[File] masks = make_mask_and_diff.mask_file
		Array[File] diffs = make_mask_and_diff.diff
		File? tax_tree = trees.taxonium_tree
		Array[File]? fastqc_reports = FastqcWF.reports
	}
}

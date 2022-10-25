version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/sort-by-name/workflows/refprep-TB.wdl" as clockwork_ref_prepWF
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/sort-by-name/tasks/map_reads.wdl" as clckwrk_map_reads
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/sort-by-name/tasks/rm_contam.wdl" as clckwrk_rm_contam
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/sort-by-name/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/handle_odd_numbers/tasks/pull_from_SRA.wdl" as sranwrp
import "https://raw.githubusercontent.com/aofarrel/enaBrowserTools-wdl/0.0.4/tasks/enaDataGet.wdl" as ena
import "https://raw.githubusercontent.com/aofarrel/mask-by-coverage/main/mask-by-coverage.wdl" as masker
import "https://raw.githubusercontent.com/aofarrel/tb_tree/add-wdl/pipelines/make_diff.wdl" as diff

workflow myco {
	input {
		Array[String] SRA_accessions
		Int min_coverage
	}

	call clockwork_ref_prepWF.ClockworkRefPrepTB

	scatter(SRA_accession in SRA_accessions) {
		call sranwrp.pull_from_SRA_directly {
			input:
				sra_accession = SRA_accession
		}
	} # output: pull_from_SRA_directly.fastqs

	# TODO: Test if this can just be nested in the above scatter instead of being its own scatter
	scatter(data in zip(SRA_accessions, pull_from_SRA_directly.fastqs)) {
		call clckwrk_map_reads.map_reads {
			input:
				unsorted_sam = true,
				sample_name = data.left,
				reads_files = data.right,
				tarball_ref_fasta_and_index = ClockworkRefPrepTB.tar_indexd_dcontm_ref,
				ref_fasta_filename = "ref.fa"
		} # output: map_reads.mapped_reads

		call masker.make_mask_file {
			input:
				sam = map_reads.mapped_reads,
				min_coverage = min_coverage
		}

		# TODO: replace with single file TSV if possible as that is much faster to localize
		call clckwrk_rm_contam.remove_contam as remove_contamination {
			input:
				bam_in = map_reads.mapped_reads,
				tarball_metadata_tsv = ClockworkRefPrepTB.tar_indexd_dcontm_ref
		} # output: remove_contamination.decontaminated_fastq_1, remove_contamination.decontaminated_fastq_2

		call clckwrk_var_call.variant_call_one_sample as varcall {
			input:
				sample_name = map_reads.mapped_reads,
				ref_dir = ClockworkRefPrepTB.tar_indexd_H37Rv_ref,
				reads_files = [remove_contamination.decontaminated_fastq_1, remove_contamination.decontaminated_fastq_2]
		} # output: varcall.vcf_final_call_set

		call diff.make_diff as diffmaker {
			input:
				vcf = varcall.vcf_final_call_set
		}
	}

	output {
		# outputting everything for debugging purposes
		Array[File] sams  = map_reads.mapped_reads
		Array[File] masks = make_mask_file.mask_file
		Array[File] dcnfq1= remove_contamination.decontaminated_fastq_1
		Array[File] dcnfq2= remove_contamination.decontaminated_fastq_2
		Array[File] minos = varcall.vcf_final_call_set
		Array[File] diffs = diffmaker.diff
	}
}
version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/1.0.0/tasks/map_reads.wdl" as clckwrk_map_reads
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/main/tasks/rm_contam.wdl" as clckwrk_rm_contam
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/main/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/main/tasks/pull_from_SRA.wdl" as sranwrp
import "https://raw.githubusercontent.com/aofarrel/mask-by-coverage/main/mask-by-coverage.wdl" as masker

workflow myco {
	input {
		File tarball_decontaminated_ref
		File tarball_H37Rv_ref
		File contamination_metadata_tsv
		Array[String] SRA_accessions
		Int min_coverage = 1
	}

	scatter(SRA_accession in SRA_accessions) {

		call sranwrp.pull_from_SRA_directly {
			input:
				sra_accession = SRA_accession

		} # output: pull_from_SRA_directly.fastqs

		call clckwrk_map_reads.map_reads {
			input:
				sample_name = SRA_accession,
				tarball_ref_fasta_and_index = tarball_decontaminated_ref,
				ref_fasta_filename = "ref.fa",
				reads_files = pull_from_SRA_directly.fastqs

		} # output: map_reads.mapped_reads

		call clckwrk_rm_contam.remove_contam {
			input:
				bam_in = map_reads.mapped_reads,
				metadata_tsv = contamination_metadata_tsv
		} # output: remove_contam.decontaminated_fastq_1, remove_contam.decontaminated_fastq_2

		call masker.make_mask_file {
			input:
				sam = map_reads.mapped_reads,
				min_coverage = min_coverage
		}

		call clckwrk_var_call.variant_call_one_sample {
			input:
				sample_name = map_reads.mapped_reads,
				ref_dir = tarball_H37Rv_ref,
				reads_files = [remove_contam.decontaminated_fastq_1, remove_contam.decontaminated_fastq_2]
		}
	}
	
}
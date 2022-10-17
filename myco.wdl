version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/sort-by-name/workflows/refprep-TB.wdl" as clockwork_ref_prepWF
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/sort-by-name/tasks/map_reads.wdl" as clckwrk_map_reads
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/sort-by-name/tasks/rm_contam.wdl" as clckwrk_rm_contam
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/sort-by-name/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/main/tasks/pull_from_SRA.wdl" as sranwrp
import "https://raw.githubusercontent.com/aofarrel/enaBrowserTools-wdl/0.0.4/tasks/enaDataGet.wdl" as ena
import "https://raw.githubusercontent.com/aofarrel/mask-by-coverage/main/mask-by-coverage.wdl" as masker

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

	scatter(data in zip(SRA_accessions, pull_from_SRA_directly.fastqs)) {
		call clckwrk_map_reads.map_reads {
			input:
				sample_name = data.left,
				reads_files = data.right,
				tarball_ref_fasta_and_index = ClockworkRefPrepTB.tar_indexd_dcontm_ref,
				ref_fasta_filename = "ref.fa",
				reads_files = pull_from_SRA_directly.fastqs
		}

	} # output: map_reads.mapped_reads

	scatter(sam_file in map_reads.mapped_reads) {
		# TODO: replace with single file TSV if possible as that is much faster to localize
		call clckwrk_rm_contam.remove_contam as remove_contamination {
			input:
				bam_in = sam_file,
				tarball_metadata_tsv = ClockworkRefPrepTB.tar_indexd_dcontm_ref,
		} # output: remove_contam.decontaminated_fastq_1, remove_contam.decontaminated_fastq_2

		call masker.make_mask_file {
			input:
				sam = map_reads.mapped_reads,
				min_coverage = min_coverage
		}

		call clckwrk_var_call.variant_call_one_sample {
			input:
				sample_name = sam_file,
				ref_dir = ClockworkRefPrepTB.tar_indexd_H37Rv_ref,
				reads_files = [remove_contamination.decontaminated_fastq_1, remove_contamination.decontaminated_fastq_2]
		}
	}
}
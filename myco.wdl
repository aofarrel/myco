version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/add-var-call-debugging-task/workflows/refprep-TB.wdl" as clockwork_ref_prepWF
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/main/tasks/map_reads.wdl" as clckwrk_map_reads
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/main/tasks/rm_contam.wdl" as clckwrk_rm_contam
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/add-var-call-debugging-task/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/main/tasks/pull_from_SRA.wdl" as sranwrp
import "https://raw.githubusercontent.com/aofarrel/enaBrowserTools-wdl/0.0.4/tasks/enaDataGet.wdl" as ena
import "https://raw.githubusercontent.com/aofarrel/mask-by-coverage/main/mask-by-coverage.wdl" as masker
import "https://raw.githubusercontent.com/aofarrel/tb_tree/add-wdl/pipelines/make_diff.wdl" as diff

workflow myco {
	input {
		Array[String] SRA_accessions
		Int min_coverage
		Boolean? skip_decontamination = false
	}

	call clockwork_ref_prepWF.ClockworkRefPrepTB

	scatter(SRA_accession in SRA_accessions) {
		call sranwrp.pull_from_SRA_directly {
			input:
				sra_accession = SRA_accession
		} # output: pull_from_SRA_directly.fastqs
		
		if(length(pull_from_SRA_directly.fastqs)>1) {
				Array[File] paired_fastqs=select_all(pull_from_SRA_directly.fastqs)
		}
	}
	
	Array[Array[File]] pulled_fastqs = select_all(paired_fastqs)

	if(!skip_decontamination) {
		scatter(pulled_fastq in pulled_fastqs) {
			call clckwrk_map_reads.map_reads as map_reads_for_decontam {
				input:
					unsorted_sam = true,
					reads_files = pulled_fastq,
					tarball_ref_fasta_and_index = ClockworkRefPrepTB.tar_indexd_dcontm_ref,
					ref_fasta_filename = "ref.fa"
			} # output: map_reads_for_decontam.mapped_reads

			# TODO: replace with single file TSV if possible as that is much faster to localize
			call clckwrk_rm_contam.remove_contam as remove_contamination {
				input:
					bam_in = map_reads_for_decontam.mapped_reads,
					tarball_metadata_tsv = ClockworkRefPrepTB.tar_indexd_dcontm_ref
			} # output: remove_contamination.decontaminated_fastq_1, remove_contamination.decontaminated_fastq_2

			call clckwrk_var_call.variant_call_one_sample_cool as varcall {
				input:
					sample_name = map_reads_for_decontam.mapped_reads,
					ref_dir = ClockworkRefPrepTB.tar_indexd_H37Rv_ref,
					reads_files = [remove_contamination.decontaminated_fastq_1, remove_contamination.decontaminated_fastq_2]
			} # output: varcall.vcf_final_call_set, varcall.mapped_to_ref

		}
	}

	if(skip_decontamination) {
		scatter(pulled_fastq in pulled_fastqs) {
			call clckwrk_var_call.variant_call_one_sample_cool as varcall_no_decontam {
				input:
					sample_name = map_reads_for_decontam.mapped_reads,
					ref_dir = ClockworkRefPrepTB.tar_indexd_H37Rv_ref,
					reads_files = pulled_fastq
			} # output: varcall.vcf_final_call_set, varcall.mapped_to_ref

		}
	}

	Array[File] minos_vcfs=select_all(select_first([varcall_no_decontam.vcf_final_call_set, varcall.vcf_final_call_set])
	Array[File] bams_to_ref=select_all(select_first([varcall_no_decontam.mapped_to_ref, varcall.mapped_to_ref)

	scatter(bam_to_ref in bams_to_ref) {
		call masker.make_mask_file {
			input:
				bam = bam_to_ref,
				min_coverage = min_coverage
		}
	}


	scatter(minos_vcf in minos_vcfs) {
		call diff.make_diff as diffmaker {
			input:
				vcf = minos_vcf
		}
	}

	output {
		# outputting everything for debugging purposes
		Array[File] reads_mapped_to_decontam  = map_reads_for_decontam.mapped_reads
		Array[File] reads_mapped_to_H37Rv = bams_to_ref
		Array[File] masks = make_mask_file.mask_file
		Array[File] dcnfq1= remove_contamination.decontaminated_fastq_1
		Array[File] dcnfq2= remove_contamination.decontaminated_fastq_2
		Array[File] minos = minos_vcfs
		Array[File] diffs = diffmaker.diff
		Array[File?] debug_error = varcall.debug_error
	}
}
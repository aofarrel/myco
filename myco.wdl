version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.0.0/workflows/refprep-TB.wdl" as clockwork_ref_prepWF
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.0.0/tasks/map_reads.wdl" as clckwrk_map_reads
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.0.0/tasks/rm_contam.wdl" as clckwrk_rm_contam
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.0.0/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/main/tasks/pull_fastqs.wdl" as sranwrp_pull
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/main/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/enaBrowserTools-wdl/0.0.4/tasks/enaDataGet.wdl" as ena
import "https://raw.githubusercontent.com/aofarrel/mask-by-coverage/main/mask-by-coverage.wdl" as masker
import "https://raw.githubusercontent.com/aofarrel/parsevcf/main/vcf_to_diff.wdl" as diff

workflow myco {
	input {
		File biosample_accessions
		File typical_tb_masked_regions
		Int min_coverage
		Boolean skip_decontamination = false
	}

	call clockwork_ref_prepWF.ClockworkRefPrepTB

	call sranwrp_processing.extract_accessions_from_file as get_sample_IDs {
		input:
			accessions_file = biosample_accessions
	}

	scatter(biosample_accession in get_sample_IDs.accessions) {
		call sranwrp_pull.pull_fq_from_biosample as pull {
			input:
				biosample_accession = biosample_accession
		} # output: pull.fastqs

		if(length(pull.fastqs)>1) {
    		Array[File] paired_fastqs=select_all(pull.fastqs)
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

			call clckwrk_var_call.variant_call_one_sample_verbose as varcall {
				input:
					sample_name = map_reads_for_decontam.mapped_reads,
					ref_dir = ClockworkRefPrepTB.tar_indexd_H37Rv_ref,
					reads_files = [remove_contamination.decontaminated_fastq_1, remove_contamination.decontaminated_fastq_2]
			} # output: varcall.vcf_final_call_set, varcall.mapped_to_ref

		}
	}

	if(skip_decontamination) {
		scatter(pulled_fastq in pulled_fastqs) {
			call clckwrk_var_call.variant_call_one_sample_verbose as varcall_no_decontam {
				input:
					sample_name = basename(pulled_fastq[0], ".fq"),
					ref_dir = ClockworkRefPrepTB.tar_indexd_H37Rv_ref,
					reads_files = pulled_fastq
			} # output: varcall.vcf_final_call_set, varcall.mapped_to_ref

		}
	}

	Array[File] minos_vcfs=select_all(select_first([varcall_no_decontam.vcf_final_call_set, varcall.vcf_final_call_set]))
	Array[File] bams_to_ref=select_all(select_first([varcall_no_decontam.mapped_to_ref, varcall.mapped_to_ref]))


	scatter(vcfs_and_bams in zip(bams_to_ref, minos_vcfs)) {
		call diff.make_mask_and_diff {
			input:
				bam = vcfs_and_bams.left,
				vcf = vcfs_and_bams.right,
				min_coverage = min_coverage,
				tbmf = typical_tb_masked_regions
		}
	}

	output {
		# outputting everything for debugging purposes
		Array[File]? reads_mapped_to_decontam  = map_reads_for_decontam.mapped_reads
		Array[File] reads_mapped_to_H37Rv = bams_to_ref
		Array[File] masks = make_mask_file.mask_file
		Array[File]? dcnfq1= remove_contamination.decontaminated_fastq_1
		Array[File]? dcnfq2= remove_contamination.decontaminated_fastq_2
		Array[File] minos = minos_vcfs
		Array[File] diffs = diffmaker.diff
		Array[File?] debug_error = select_first([varcall.debug_error, varcall_no_decontam.debug_error])
	}
}
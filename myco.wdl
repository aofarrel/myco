version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.0.1/workflows/refprep-TB.wdl" as clockwork_ref_prepWF
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.3.1/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.3.1/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.3/tasks/pull_fastqs.wdl" as sranwrp_pull
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.2/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/parsevcf/0.99.1/vcf_to_diff.wdl" as diff

workflow myco {
	input {
		File biosample_accessions
		File typical_tb_masked_regions
		Int min_coverage = 10
		Boolean tar_fqs = false
		Boolean less_scattering = false
	}

	call clockwork_ref_prepWF.ClockworkRefPrepTB

	call sranwrp_processing.extract_accessions_from_file as get_sample_IDs {
		input:
			accessions_file = biosample_accessions
	}

	scatter(biosample_accession in get_sample_IDs.accessions) {
		call sranwrp_pull.pull_fq_from_biosample as pull {
			input:
				biosample_accession = biosample_accession,
				tar_outputs = tar_fqs
		} # output: pull.fastqs OR pull.tarball_fastqs
		if(length(pull.fastqs)>1) {
    		Array[File] paired_fastqs=select_all(pull.fastqs)
  		}
	}

	call sranwrp_processing.cat_files as cat_reports {
		input:
			files = pull.results
	}

	if(!less_scattering) {
		Array[Array[File]] pulled_fastqs   = select_all(paired_fastqs)
		scatter(pulled_fastq in pulled_fastqs) {
			call clckwrk_combonation.combined_decontamination_single as decontaminate_one_sample {
				input:
					unsorted_sam = true,
					reads_files = pulled_fastq,
					tarball_ref_fasta_and_index = ClockworkRefPrepTB.tar_indexd_dcontm_ref,
					ref_fasta_filename = "ref.fa"
			}

			call clckwrk_var_call.variant_call_one_sample_simple as varcall_with_array {
				input:
					ref_dir = ClockworkRefPrepTB.tar_indexd_H37Rv_ref,
					reads_files = [decontaminate_one_sample.decontaminated_fastq_1, decontaminate_one_sample.decontaminated_fastq_2]
			} # output: varcall_with_array.vcf_final_call_set, varcall_with_array.mapped_to_ref

		}

		Array[File] minos_vcfs_=select_all(varcall_with_array.vcf_final_call_set)
		Array[File] bams_to_ref_=select_all(varcall_with_array.mapped_to_ref)


		scatter(vcfs_and_bams in zip(bams_to_ref_, minos_vcfs_)) {
			call diff.make_mask_and_diff as make_mask_and_diff_ {
				input:
					bam = vcfs_and_bams.left,
					vcf = vcfs_and_bams.right,
					min_coverage = min_coverage,
					tbmf = typical_tb_masked_regions
			}
		}
	}

	if(less_scattering) {
		Array[File] tarball_paired_fastqs=select_all(pull.tarball_fastqs)
		call clckwrk_combonation.combined_decontamination_multiple as decontaminate_many_samples {
			input:
				unsorted_sam = true,
				tarballs_of_read_files = tarball_paired_fastqs,
				tarball_ref_fasta_and_index = ClockworkRefPrepTB.tar_indexd_dcontm_ref,
				ref_fasta_filename = "ref.fa"
		} # output: decontaminate_many_samples.tarballs_of_decontaminated_reads
		
		scatter(one_sample in decontaminate_many_samples.tarballs_of_decontaminated_reads) {
			call clckwrk_var_call.variant_call_one_sample_verbose as varcall_with_tarballs {
				input:
					ref_dir = ClockworkRefPrepTB.tar_indexd_H37Rv_ref,
					tarball_of_reads_files = one_sample
			} # output: varcall_with_tarballs.vcf_final_call_set, varcall_with_tarballs.mapped_to_ref
		}

		Array[File] minos_vcfs=select_all(varcall_with_tarballs.vcf_final_call_set)
		Array[File] bams_to_ref=select_all(varcall_with_tarballs.mapped_to_ref)


		scatter(vcfs_and_bams in zip(bams_to_ref, minos_vcfs)) {
			call diff.make_mask_and_diff as make_mask_and_diff {
				input:
					bam = vcfs_and_bams.left,
					vcf = vcfs_and_bams.right,
					min_coverage = min_coverage,
					tbmf = typical_tb_masked_regions
			}
		}

	}

	output {
		Array[File] minos = select_first([minos_vcfs, minos_vcfs_])
		Array[File] masks = select_first([make_mask_and_diff.mask_file, make_mask_and_diff_.mask_file])
		Array[File] diffs = select_first([make_mask_and_diff.diff, make_mask_and_diff_.diff])
		File pull_report = cat_reports.outfile
	}
}

version 1.0
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.16.5/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/fastp-wdl/0.0.4/fastp_tasks.wdl" as fastp
import "https://raw.githubusercontent.com/aofarrel/vcf_to_diff_wdl/0.0.3/vcf_to_diff.wdl" as diff

# This is a stripped-down version of myco which only runs the variant caller, vcf-to-diff, and (optionally) fastp.
# There is NO decontamination, TB-Profiler, covstats, or QC beyond "does the variant caller crash/time out or not."
# Need to decontaminate your fastqs, but don't want to run myco_raw? Try Decontam_And_Combine_One_Samples_Fastqs, which can be found
# on Dockstore: https://dockstore.org/workflows/github.com/aofarrel/clockwork-wdl/Decontam_And_Combine_One_Samples_Fastqs

/*
 * Copyright (C) 2025 Ash O'Farrell
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */



workflow myco {
	input {
		Array[Array[File]] paired_decontaminated_fastq_sets
		
		Boolean fastp_clean            = true
		Float   max_low_coverage_sites =   0.05
		Int     min_coverage_per_site  =  10
		Int     timeout_variant_caller = 120
	}

	parameter_meta {
		fastp_clean: "If true, clean reads with fastp before calling variants"
		max_low_coverage_sites: "If a diff file has higher than this percent (0.5 = 50%) bad data, do not include it in the tree"
		min_coverage_per_site: "Positions with coverage below this value will be masked in diff files"
		paired_decontaminated_fastq_sets: "Nested array of decontaminated and merged fastq pairs. Each inner array represents one sample; each sample needs precisely one gzipped forward read and one gzipped reverse read."
		timeout_variant_caller: "Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to 0 to never timeout)"
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
			}
	}

	output {
		Array[File] vcfs = minos_vcfs
		Array[File] masks = make_mask_and_diff.mask_file
		Array[File?] diffs = make_mask_and_diff.diff
	}
}

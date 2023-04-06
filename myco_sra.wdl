version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.8.0/workflows/refprep-TB.wdl" as clockwork_ref_prepWF
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.8.0/tasks/combined_decontamination.wdl" as clckwrk_combonation
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/2.8.0/tasks/variant_call_one_sample.wdl" as clckwrk_var_call
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.8/tasks/pull_fastqs.wdl" as sranwrp_pull
import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.8/tasks/processing_tasks.wdl" as sranwrp_processing
import "https://raw.githubusercontent.com/aofarrel/tree_nine/0.0.5/tree_nine.wdl" as build_treesWF
import "https://raw.githubusercontent.com/aofarrel/parsevcf/1.1.4/vcf_to_diff.wdl" as diff
import "https://raw.githubusercontent.com/aofarrel/fastqc-wdl/main/fastqc.wdl" as fastqc
import "https://raw.githubusercontent.com/aofarrel/tb_profiler/0.1.1/tb_profiler.wdl" as profiler

workflow myco {
	input {
		File biosample_accessions
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
		Int     timeout_decontam_part1 =   20
		Int     timeout_decontam_part2 =   15
		Int     timeout_variant_caller =  120
	}

	parameter_meta {
		biosample_accessions: "File of BioSample accessions to pull, one accession per line"
		bad_data_threshold: "If a diff file has higher than this percent (0.5 = 50%) bad data, do not include it in the tree"
		decorate_tree: "Should usher, taxonium, and NextStrain trees be generated? Requires input_tree and ref_genome"
		fastqc_on_timeout: "If true, fastqc one read from a sample when decontamination or variant calling times out"
		force_diff: "If true and if decorate_tree is false, generate diff files. (Diff files will always be created if decorate_tree is true.)"
		input_tree: "Base tree to use if decorate_tree = true"
		min_coverage: "Positions with coverage below this value will be masked in diff files"
		ref_genome_for_tree_building: "Ref genome for building trees iff different from ref genome used to call variants -- must have ONLY `>NC_000962.3` on its first line"
		subsample_cutoff: "If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)"
		subsample_seed: "Seed used for subsampling with seqtk"
		timeout_decontam_part1: "Discard any sample that is still running in clockwork map_reads after this many minutes (set to 0 to never timeout)"
		timeout_decontam_part2: "Discard any sample that is still running in clockwork rm_contam after this many minutes (set to 0 to never timeout)"
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

	call sranwrp_processing.extract_accessions_from_file as get_sample_IDs {
		input:
			accessions_file = biosample_accessions,
			filter_na = true
	}

	scatter(biosample_accession in get_sample_IDs.accessions) {
		call sranwrp_pull.pull_fq_from_biosample as pull {
			input:
				biosample_accession = biosample_accession,
				fail_on_invalid = false,
				subsample_cutoff = subsample_cutoff,
				subsample_seed = subsample_seed,
				tar_outputs = false
		} # output: pull.fastqs
		if(length(pull.fastqs)>1) {
    		Array[File] paired_fastqs=select_all(pull.fastqs)
  		}
	}

	call sranwrp_processing.cat_strings as cat_reports {
		input:
			strings = pull.results
	}

	Array[Array[File]] pulled_fastqs = select_all(paired_fastqs)
	scatter(pulled_fastq in pulled_fastqs) {
		call clckwrk_combonation.combined_decontamination_single as per_sample_decontam {
			input:
				unsorted_sam = true,
				reads_files = pulled_fastq,
				tarball_ref_fasta_and_index = ClockworkRefPrepTB.tar_indexd_dcontm_ref,
				ref_fasta_filename = "ref.fa",
				filename_metadata_tsv = "remove_contam_metadata.tsv",
				timeout_map_reads = timeout_decontam_part1,
				timeout_decontam = timeout_decontam_part2
		}

		if(defined(per_sample_decontam.decontaminated_fastq_1)) {
		# This region only executes if decontaminated fastqs exist.
		# We can use this to coerce File? into File by using a
		# select_first() where the first element is the File? we know
		# absolutely must exist, and the second element is bogus
    		File real_decontaminated_fastq_1=select_first([
    			per_sample_decontam.decontaminated_fastq_1, 
    				biosample_accessions])
    		File real_decontaminated_fastq_2=select_first(
    			[per_sample_decontam.decontaminated_fastq_2, 
    				biosample_accessions])

			call clckwrk_var_call.variant_call_one_sample_simple as varcall_with_array {
				input:
					ref_dir = ClockworkRefPrepTB.tar_indexd_H37Rv_ref,
					reads_files = [real_decontaminated_fastq_1, real_decontaminated_fastq_2],
					timeout = timeout_variant_caller
			}
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

		call profiler.tb_profiler_bam as profile {
			input:
				bam = vcfs_and_bams.left
		}
	}

	if(fastqc_on_timeout) {
		Array[File] bad_fastqs_decontam_ = select_all(per_sample_decontam.check_this_fastq)
		Array[File] bad_fastqs_varcallr_ = select_all(varcall_with_array.check_this_fastq)
		Array[Array[File]] bad_fastqs_   = [bad_fastqs_decontam_, bad_fastqs_varcallr_]
		if(length(per_sample_decontam.check_this_fastq)>1 && length(bad_fastqs_varcallr_)>1) {
			Array[File] bad_fastqs_both      = flatten(bad_fastqs_)  
		}
		if(length(per_sample_decontam.check_this_fastq)>1) {
			Array[File] bad_fastqs_decontam = select_all(per_sample_decontam.check_this_fastq)
		}
		if(length(bad_fastqs_varcallr_)>1) {
			Array[File] bad_fastqs_varcallr = select_all(bad_fastqs_varcallr_)
		}
		call fastqc.FastqcWF {
			input:
				fastqs = select_first([bad_fastqs_both, bad_fastqs_decontam, bad_fastqs_varcallr])
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
				ref = select_first([ref_genome_for_tree_building, ClockworkRefPrepTB.reference_genome_fasta]),
				coverage_reports = coerced_reports,
				bad_data_threshold = bad_data_threshold
		}
	}

	output {
		File download_report = cat_reports.outfile
		Array[File] minos = minos_vcfs
		Array[File] masks = make_mask_and_diff.mask_file
		Array[File] tbprofiler = profile.tbprofiler_results
		Array[File?] diffs = make_mask_and_diff.diff
		File? tree_usher = trees.usher_tree
		File? tree_taxonium = trees.taxonium_tree
		File? tree_nextstrain = trees.nextstrain_tree
		Array[File]? trees_nextstrain = trees.nextstrain_subtrees
		Array[File]? fastqc_reports = FastqcWF.reports
	}
}

version 1.0

import "./myco_cleaned.wdl" as WF

# This is just a one-sample wrapper for myco_cleaned. It is intended for Terra data tables with a format like this:







workflow myco_cleaned_one_sample {
    input {
        File decontaminated_fastq_1
        File decontaminated_fastq_2
    }

    call WF.myco {
        input:
            paired_fastq_sets = [[decontaminated_fastq_1, decontaminated_fastq_2]]
    }

output {
		Array[File] minos = WF.myco.minos_vcfs
		Array[File] masks = WF.myco.make_mask_and_diff.mask_file
		Array[File?] diffs = WF.myco.make_mask_and_diff.diff
		File? tax_tree = WF.myco.trees.taxonium_tree
		Array[File]? fastqc_reports = WF.myco.FastqcWF.reports
	}
}

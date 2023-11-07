version 1.0

# Replace the first URL here with the URL of the workflow to be checked.
import "https://raw.githubusercontent.com/aofarrel/checker-WDL-templates/v1.1.0/check_wf_outputs/outputs_some_optional/parent_opt.wdl" as myco_raw
import "https://raw.githubusercontent.com/dockstore/checker-WDL-templates/v1.1.0/checker_tasks/filecheck_task.wdl" as verify_file
import "https://raw.githubusercontent.com/dockstore/checker-WDL-templates/v1.1.0/checker_tasks/arraycheck_task.wdl" as verify_array
import "./myco_raw.wdl" as myco_rawest # "myco_raw" causes a namespace conflict

workflow checker {
	input {
		Array[Array[File]] paired_fastq_sets

		
		Array[File] truth_myco_raw_default # vcf, bam, mask, diff, qc csv
		File truth_myco_raw_strict # just the qc csv
	}

	# Default settings
	call myco_rawest.myco as myco_raw_default {
		input:
			paired_fastq_sets = paired_fastq_sets
	}
	
	call verify_array.arraycheck_classic as check_myco_raw_default {
		input:
			test = [myco_raw_default.vcfs[0], myco_raw_default.bams[0], myco_raw_default.masks[0], myco_raw_default.diffs[0], myco_raw_default.qc_csv],
			truth = arrayTruth
	}
	
    # Extremely strict QC
    call myco_rawest.myco as myco_raw_strict {
		input:
			paired_fastq_sets = paired_fastq_sets
			# TODO: better define values for this!
	}

	call verify_file.filecheck as check_myco_raw_strict {
		input:
			test = myco_raw_strict.qc_csv,
			truth = truth_myco_raw_strict
	}
	
	# Extremely strict QC (basically a forced failure)
    call myco_rawest.myco as myco_raw_super_strict {
		input:
			paired_fastq_sets = paired_fastq_sets,
			tbprofilerQC_min_pct_mapped = 100,
			soft_pct_mapped = false
	}

	call verify_file.filecheck as check_myco_raw_super_strict {
		input:
			test = myco_raw_super_strict.qc_csv,
			truth = truth_myco_raw_strict
	}

}
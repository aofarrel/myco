version 1.0
import "https://raw.githubusercontent.com/dockstore/checker-WDL-templates/v1.1.0/checker_tasks/arraycheck_task.wdl" as verify_array
import "./myco_raw.wdl" as myco_raw

# Unlike myco raw, this checker workflow is designed to run as seperate workflow instances.
# To restate: To run myco_raw on n samples, you can either:
# a) Use a Terra data table where each row represents one sample --> launches n instances of 
#    myco_raw where each instance is only 'aware' of a single sample
# b) Define multiple inner arrays where each inner array is one sample's fastqs --> launches 
#    one instance of myco_raw and that instance is 'aware' of all n samples
# This checker workflow IS NOT LIKE THAT! Only option (a) will work as expected. A Terra data
# table with open-access data and truth files has been provided for this reason.

workflow checker {
	input {
		# inner array is one sample's fastqs. length of inner array can be any even number above 0.
		# length of outer array can only be one (since this is a single-sample workflow).
		Array[Array[File]] paired_fastq_sets
		
		String TRUTH_code
		Array[File] FALLBACK_bai
		Array[File] FALLBACK_diff
		Array[File] FALLBACK_report
		Array[File] FALLBACK_vcf
		
		# These are arrays, but they should only contain one (or zero) values
		Array[File] TRUTH_bai
		Array[File] TRUTH_diff
		Array[File] TRUTH_qc
		Array[File] TRUTH_report  # TODO fix or remove report
		Array[File] TRUTH_vcf
	}

	# Default settings -- this should also catch if any meaningful default settings changed
	call myco_raw.myco as myco_raw_default {
		input:
			paired_fastq_sets = paired_fastq_sets
	}
	
	Array[File] TEST_bai = flatten(select_all([myco_raw_default.bais, TRUTH_bai, FALLBACK_bai]))
	Array[File] TEST_diff = flatten(select_all([myco_raw_default.diffs, TRUTH_diff, FALLBACK_diff]))
	File TEST_qc = myco_raw_default.qc_csv
	Array[File] TEST_report = select_all(select_first([myco_raw_default.diff_reports, TRUTH_report, FALLBACK_report])) # awkward due to being Array[File?]
	Array[File] TEST_vcf = flatten(select_all([myco_raw_default.vcfs, TRUTH_vcf, FALLBACK_vcf]))
	
	call verify_array.arraycheck_classic as check_myco_raw_default {
		input:
			test = flatten([TEST_bai, TEST_diff, [TEST_qc], TEST_report, TEST_vcf]),
			truth = flatten(select_all([TRUTH_bai, TRUTH_diff, TRUTH_qc, TRUTH_report, TRUTH_vcf]))
	}
	
    # Extremely strict QC
    #call myco_raw.myco as myco_raw_strict {
	#	input:
	#		paired_fastq_sets = paired_fastq_sets
	#		# TODO: define what strict actually means!
	#}

	#call verify_file.filecheck as check_myco_raw_strict {
	#	input:
	#		test = myco_raw_strict.qc_csv,
	#		truth = truth_myco_raw_strict
	#}
	
	# Extremely strict QC (basically a forced failure)
    #call myco_raw.myco as myco_raw_super_strict {
	#	input:
	#		paired_fastq_sets = paired_fastq_sets,
	#		tbprofilerQC_min_pct_mapped = 100,
	#		soft_pct_mapped = false
	#}

	#call verify_file.filecheck as check_myco_raw_super_strict {
	#	input:
	#		test = myco_raw_super_strict.qc_csv,
	#		truth = truth_myco_raw_strict
	#}

}
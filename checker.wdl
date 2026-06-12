version 1.0
import "https://raw.githubusercontent.com/aofarrel/checker-WDL-templates/disk-size-override/checker_tasks/arraycheck_task.wdl" as verify_array
import "https://raw.githubusercontent.com/aofarrel/checker-WDL-templates/add-stringcheck/checker_tasks/stringcheck_task.wdl" as verify_string
import "./myco_raw.wdl" as myco_raw
import "./myco_sra.wdl" as myco_sra

# This workflow is designed to run on a data table in Terra. You can run it locally, but you'll need
# to run it each row of checker_data_table as a separate workflow.
#
# Terra instructions:
# 1) Import this checker workflow
# 2) Import the checker_data_table as a TSV to create a Terra data table
# 3) In workflows tab, select checker workflow
# 4) Select "Run workflow(s) with inputs defined by data table"
# 5) Select the checker workflow data table
# 6) Fill in inputs
# 7) Optional: Set a workflow cost threshold
# 8) Launch

workflow checker {
	input {
		
		# In Terra: myco_raw_paired_fastq_sets = [this.myco_raw_input_fqs]
		Array[Array[File]] myco_raw_paired_fastq_sets

		# In Terra: myco_sra_biosample = this.myco_sra_input_BioSample
		String myco_sra_biosample

		# Fallback file for known QC-failing samples. A QC-failing sample will give no output,
		# so in order to run a checker on it, select_first([expected_output, fallback]). In the
		# expected success case, the comparison will be against expected_output and the actual
		# pipeline out. In the expected failure case, the comparison will be against two
		# instances of the fallback file. If it fails when success expected or succeeds when
		# fail expected, one will be the fallback file and the other won't, resulting in (correct)
		# error to flag the mismatch.
		File fallback
		
		# These are arrays, but except for myco_raw multi sample test cases (currently not included),
		# they should only contain one value. If expecting no file (ie known QC fail) then that file 
		# will be the fallback file.

		# When running myco_raw at default values, these should be the result
		Array[File] TRUTH_mycoraw_default_diff
		Array[File] TRUTH_mycoraw_default_diff_report
		Array[File] TRUTH_mycoraw_default_decontam_report
		String TRUTH_mycoraw_default_status

		# When running myco_sra at default values, these should be the result
		Array[File] TRUTH_mycosra_default_diff
		Array[File] TRUTH_mycosra_default_diff_report
		Array[File] TRUTH_mycosra_default_decontam_report
		String TRUTH_mycosra_default_status

		# When running myco_raw at default values, except just_like_2024 = true
		# TODO: Compare these outputs to the actual "publication reproducibility" branch
		#Array[File] TRUTH_mycoraw_legacy_clockworkCHM13decon_bai
		#Array[File] TRUTH_mycoraw_legacy_clockworkCHM13decon_diff
		#Array[File] TRUTH_mycosra_legacy_clockworkCHM13decon_bai
		#Array[File] TRUTH_mycosra_legacy_clockworkCHM13decon_diff

		# When running myco_raw at default values, except just_like_2024 = true and use_varpipe = true
		# Deprioritized, nobody uses varpipe
		#Array[File] TRUTH_legacy_varpipedecon_bai
		#Array[File] TRUTH_legacy_varpipedecon_diff

		Int checker_disk_size_override
	}

	# Turn fallback into a one-element array so it can be used as a fallback for constructing Array[File]
	Array[File] fallback_array = [fallback]

	# Default settings -- this should also catch if any meaningful default settings changed
	call myco_raw.myco as myco_raw_default {
		input:
			paired_fastq_sets = myco_raw_paired_fastq_sets
	}

	# Most of these are considered Array[File] (even though they can be empty) so we can just select_first()
	# tbd_decontam_reports is considered Array[File?] so it needs to be chained with select_all first
	
	Array[File] TEST_mycoraw_default_diff = select_first([myco_raw_default.tbd_diffs, fallback_array])
	Array[File] TEST_mycoraw_default_diff_report = select_first([select_all(myco_raw_default.tbd_diff_reports), fallback_array])
	Array[File] TEST_mycoraw_default_decontam_report = select_first([select_all(myco_raw_default.tbd_decontam_reports), fallback_array])

	call verify_string.stringcheck as status_myco_raw_default {
		input:
			test = myco_raw_default.tbd_status,
			truth = TRUTH_mycoraw_default_status
	}

	call verify_array.arraycheck_classic as diff_myco_raw_default {
		input:
			test = TEST_mycoraw_default_diff,
			truth = TRUTH_mycoraw_default_diff,
			disk_size_override = checker_disk_size_override
	}

	call verify_array.arraycheck_classic as diffreport_myco_raw_default {
		input:
			test = TEST_mycoraw_default_diff_report,
			truth = TRUTH_mycoraw_default_diff_report,
			disk_size_override = checker_disk_size_override
	}

	call verify_array.arraycheck_classic as decontam_myco_raw_default {
		input:
			test = TEST_mycoraw_default_decontam_report,
			truth = TRUTH_mycoraw_default_decontam_report,
			disk_size_override = checker_disk_size_override
	}

	# Now do myco_sra
	call myco_sra.myco as myco_sra_default {
		input:
			biosample_accession_str = myco_sra_biosample
	}

	Array[File] TEST_mycosra_default_diff = select_first([myco_sra_default.tbd_diffs, fallback_array])
	Array[File] TEST_mycosra_default_diff_report = select_first([select_all(myco_sra_default.tbd_diff_reports), fallback_array])
	Array[File] TEST_mycosra_default_decontam_report = select_first([select_all(myco_sra_default.tbd_decontam_reports), fallback_array])

	call verify_string.stringcheck as status_myco_sra_default {
		input:
			test = myco_sra_default.tbd_status,
			truth = TRUTH_mycosra_default_status
	}

	call verify_array.arraycheck_classic as diff_myco_sra_default {
		input:
			test = TEST_mycosra_default_diff,
			truth = TRUTH_mycosra_default_diff,
			disk_size_override = checker_disk_size_override
	}

	call verify_array.arraycheck_classic as diffreport_myco_sra_default {
		input:
			test = TEST_mycosra_default_diff_report,
			truth = TRUTH_mycosra_default_diff_report,
			disk_size_override = checker_disk_size_override
	}

	call verify_array.arraycheck_classic as decontam_myco_sra_default {
		input:
			test = TEST_mycosra_default_decontam_report,
			truth = TRUTH_mycosra_default_decontam_report,
			disk_size_override = checker_disk_size_override
	}


	
    # Extremely strict QC
    #call myco_raw.myco as myco_raw_strict {
	#	input:
	#		paired_fastq_sets = paired_fastq_sets
	#		# TODO: define what strict actually means!
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

version 1.0

import "https://raw.githubusercontent.com/aofarrel/myco/update-myco-cleaned/myco_simple.wdl" as WF

# This is just a one-sample wrapper for myco_simple. It is intended for Terra data tables with a format like this:
#
#  | entity:sample_id | FASTQ_forward             | FASTQ_reverse             | decontaminated_fastq_1 | decontaminated_fastq_2 |
#  |------------------|---------------------------|---------------------------|------------------------|------------------------|
#  | sampleA          | raw_A_r1.fq               | raw_A_r2.fq               | decontam_A_r1.fq       | decontam_A_r2.fq       |
#  | sampleB          | raw_B1_r1.fq, rawB2_r1.fq | raw_B1_r2.fq, rawB2_r2.fq | decontam_B_r1.fq       | decontam_B_r2.fq       |
#  | sampleC          | raw_C_r1.fq               | raw_C_r2.fq               | decontam_C_r1.fq       | decontam_C_r2.fq       |
#
# However, you don't need to use this wrapper at all. You could still use the above Terra data table with myco_simple 
# 

workflow myco_simple_one_sample {
    input {
        File decontaminated_fastq_1
        File decontaminated_fastq_2
        File? typical_tb_masked_regions
    }

    call WF.myco {
        input:
            paired_decontaminated_fastq_sets = [[decontaminated_fastq_1, decontaminated_fastq_2]],
            typical_tb_masked_regions = typical_tb_masked_regions
    }

output {
		Array[File] minos = myco.minos
		Array[File] masks = myco.masks
		Array[File?] diffs = myco.diffs
		File? tree_usher = myco.tree_usher
		File? tree_taxonium = myco.tree_taxonium
		File? tree_nextstrain = myco.tree_nextstrain
		Array[File]? trees_nextstrain = myco.trees_nextstrain
		Array[File]? fastqc_reports = myco.fastqc_reports
	}
}

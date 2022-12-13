# myco: An overview

## [1] clockwork Reference Prepare
Runs my implementation of [clockwork's reference preparation standards](https://github.com/iqbal-lab-org/clockwork/wiki/Walkthrough-scripts-only#get-and-index-reference-genomes):
1. Download TB reference files
2. Index the decontamination reference
3. Index the H37Rv reference

The outputs of this process are two tarballs, with the following structure:

```
 Ref.remove_contam.tar
  ├── ref.fa
  ├── ref.fa.fai
  ├── ref.fa.minimap2_idx
  └── remove_contam_metadata.tsv

 Ref.H37Rv.tar
  ├── ref.fa
  ├── ref.fa.fai
  ├── ref.fa.minimap2_idx
  └── ref.k31.ctx
```

This is a deterministic subworkflow, and Cromwell allows for cacheing of previous workflow outputs, so this process usually only runs once. If you are using a backend that doesn't support call cacheing, you can skip this process by inputting the following:
* ClockworkRefPrepTB.bluepeter__tar_tb_ref_raw
* ClockworkRefPrepTB.bluepeter__tar_indexd_dcontm_ref
* ClockworkRefPrepTB.bluepeter__tar_indexd_H37Rv_ref

## [2] Extract BioSample accessions from input file
The user is expected to input a text containing BioSample accessions. This task grabs all unique lines in that file and outputs an Array[String] of BioSample accessions.

## scatter: each instance gets one BioSample accession

### [3] Pull fastqs for the BioSample accession
This task pulls all fastqs for a given BioSample accession using [sra-tools](https://github.com/ncbi/sra-tools). One sample might have multiple accessions; all of them are pulled. Once pulled, my script attempts to remove everything that is not a set of paired fastqs. 

For example, let's say this task got SAMN08436121. This has only one associated with it: SRR6650260. Pulling that yields three files: SRR6650260_1.fastq, SRR6650260_2.fastq, and SRR6650260.fastq. Only SRR6650260_1.fastq and SRR6650260_2.fastq will be returned.

## leave the previous scatter
There are some samples that return no valid fastqs. Because of that, we need to leave this scatter and use WDL built-in `select_all()` on our pull task's gathered output before continuing.

## [4] Decontaminate
Based on [clockwork's decontamination process](https://github.com/iqbal-lab-org/clockwork/wiki/Walkthrough-scripts-only#decontaminate-the-reads), which runs clockwork map_reads and clockwork remove_contam. I have combined these two calls into one WDL task that can run on multiple samples at once.

Running on more than one sample at a time may not always be the most economical tradeoff, so testing is currently underway to determine if it is better to scatter this task, and if so, to what extent.

## scatter: each instance gets a tarball representing one BioSample accession's decontaminated fastqs

### [5] Call variants
Based on clockwork variant_call_single, which itself combines samtools, cotex, and minos.

## leave the previous scatter

## scatter: each instance gets one bam and one vcf from task 5

### [6] Mask the outputs and create diff files
When feeding outputs into UShER, we want to make use of diff files. But first, a little bit of data processing -- it common for some regions of the TB genome to be masked. We want to avoid those problematic regions in our final output, as well as any regions without much coverage. This task cleans up our outputs and creates a diff file which can be used to make some happy little trees.

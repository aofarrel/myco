# myco: An overview

For information on running myco, [see here](./running_myco.md).

## Flowchart of myco
![Flowchart of myco](./myco_basic_flowchart.png)

## Flowchart of myco_sra
![Flowchart of myco_sra, part 1](./myco_sra_flowchart_1.png)
![Flowchart of myco_sra, part 2](./myco_sra_flowchart_2.png)

## Specific task information

### [1a] clockwork Reference Prepare
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

### [1b] Extract BioSample accessions from input file *(myco_sra.only)*
The user is expected to input a text containing BioSample accessions. This task grabs all unique lines in that file and outputs an Array[String] of BioSample accessions.

### [2] Pull fastqs for the BioSample accession *(myco_sra.only)*
This task pulls all fastqs for a given BioSample accession using [sra-tools](https://github.com/ncbi/sra-tools). One sample might have multiple accessions; all of them are pulled. Once pulled, my script attempts to remove everything that is not a set of paired fastqs. 

For example, let's say this task got SAMN08436121. This has only one run associated with it: SRR6650260. Pulling that yields three files: SRR6650260_1.fastq, SRR6650260_2.fastq, and SRR6650260.fastq. Only SRR6650260_1.fastq and SRR6650260_2.fastq will be returned.

### [3] Prepare a report of the results of the pull task *(myco_sra.only)*
There are some samples that return no valid fastqs. This task keep track of every sample's run accessions, and the result of trying to pull fastqs from each run accession.

### [4] Decontaminate
Based on [clockwork's decontamination process](https://github.com/iqbal-lab-org/clockwork/wiki/Walkthrough-scripts-only#decontaminate-the-reads), which runs clockwork map_reads and clockwork remove_contam. I have combined these two calls into one WDL task that can run on multiple samples at once. The output is a group of decontaminated fastq files.

### [5a] Call variants
Based on clockwork variant_call_single, which itself combines samtools, cotex, and minos. For each sample, the output is a single VCF file and a BAM file.

### [5b] (optional) Run FastQC on slow samples
If a sample times out in the decontamination step, it is usually due to an issue with the inputs. FastQC examines all inputs that timed out in the decontamination step, so you can see what might be going on.

### [6] Mask the outputs and create diff files
When feeding outputs into UShER, we want to make use of diff files. But first, we perform a little bit of data processing -- it common for some regions of the TB genome to be masked. We want to avoid those problematic regions in our final output, as well as any regions without much coverage. This task cleans up our outputs and creates a diff file which can be used to make some happy little trees.

### [7] (optional) Generate UShER, Taxonium, and NextStrain trees
If decorate_trees = true, and an input tree is passed in, each sample will be placed on the tree by UShER. The resulting tree will then be converted to Taxonium format, allowing it to be viewed in taxonium. NextStrain subtree JSONs will also be generated.

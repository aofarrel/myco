## Table of Contents  
    * [Workflow-level inputs](#workflow-level-inputs)
      * [FASTQ-related inputs](#fastq-related-inputs)
      * [More info on each version of myco's use case](#more-info-on-each-version-of-mycos-use-case)
      * [Non-FASTQ workflow-level inputs](#non-fastq-workflow-level-inputs)
    * [Task-level inputs](#task-level-inputs)
      * [Software settings](#software-settings)
      * [Runtime attributes](#runtime-attributes)
  
See /inputs/example_inputs.json for examples.  
  
## Workflow-level inputs  
  
### FASTQ-related inputs  
Each version of myco has a slightly different way of inputting FASTQs. A basic explanation for each workflow is in the table below. You can find more detailed explanations in each workflow's workflow-level readme.  
  
| name | type | workflow | description |  
|:---:|:---:|:---:|:---:|  
| biosample_accessions | File | myco_sra | File of BioSample accessions to pull, one accession per line |  
  
Regardless of which version of myco you use, please make sure your FASTQs:
* is Illumina paired-end data <sup>†</sup>  
* is grouped per-sample   
* len(quality scores) = len(nucleotides) for every line <sup>†</sup>  
* is actually [MTBC](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=77643)  
<sup>†</sup> myco_sra.wdl is able to detect these issues and will throw out those samples without erroring. Other forms of myco are not able to detect these issues.
It is recommend that you also keep an eye on the total size of your FASTQs. Individual files over subsample_cutoff (default450 MB, -1 disables this check) will be downsampled, but keep an eye on the cumulative size of samples. For example, a sample like SAMEA968096 has 12 run accessions associated with it. Individually, none of these run accessions' FASTQs are over 1 GB in size, but the sum total of these FASTQs could quickly fill up your disk space. (You probably should not be using SAMEA968096 anyway because it is in sample group, which can cause other issues.)

myco_cleaned expects that the FASTQs you are putting into have already been cleaned and merged. It's recommend you do this by running [Decontam_and_Combine](https://dockstore.org/workflows/github.com/aofarrel/clockwork-wdl/Decontam_And_Combine_One_Samples_Fastqs).  
  
### More info on each version of myco's use case  
* pairs of FASTQs which have been decontaminated and merged such that each sample has precisely two FASTQs associated with it**myco_cleaned** 
  * if these are in Terra data table format, you may want to use **myco_cleaned_1samp** 
 * pairs of FASTQs which have yet to be decontaminated or merged
 * if each sample has its FASTQs in a single array**myco_raw** 
 * if each sample has its forward FASTQs in one array and reverse FASTQs in another array[Decontam_And_Combine_One_Samples_Fastqs](https://dockstore.org/workflows/github.com/aofarrel/clockwork-wdl/Decontam_And_Combine_One_Samples_Fastqs), then **myco_cleaned** or **myco_cleaned_1samp** 
 * a list of SRA BioSamples whose FASTQs you'd like to use**myco_sra** 
 * a list of SRA run accessions (ERR, SRR, DRR) whose FASTQs you'd like to use[convert them to BioSamples](https://dockstore.org/workflows/github.com/aofarrel/SRANWRP/get_biosample_accessions_from_run_accessions:main?tab=info), then **myco_sra**)   
  
### Non-FASTQ workflow-level inputs  
  
| name | type | default | description |  
|:---:|:---:|:---:|:---:|  
| bad_data_threshold | Float  | 0.05 | If a diff file has higher than this percent (0.5 = 50%) bad data, do not include it in the tree |  
| decorate_tree | Boolean  | false | Should usher, taxonium, and NextStrain trees be generated? Requires input_tree and ref_genome |  
| fastqc_on_timeout | Boolean  | false | If true, fastqc one read from a sample when decontamination or variant calling times out |  
| force_diff | Boolean  | false | If true and if decorate_tree is false, generate diff files. (Diff files will always be created if decorate_tree is true.) |  
| input_tree | File? |  | Base tree to use if decorate_tree = true |  
| min_coverage | Int  | 10 | Positions with coverage below this value will be masked in diff files |  
| ref_genome_for_tree_building | File? |  | Ref genome for building trees -- must have ONLY '>NC_000962.3' on its first line |  
| subsample_cutoff | Int  | 450 | If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable) |  
| subsample_seed | Int  | 1965 | Seed used for subsampling with seqtk |  
| timeout_decontam_part1 | Int  | 20 | Discard any sample that is still running in clockwork map_reads after this many minutes (set to 0 to never timeout |  
| timeout_decontam_part2 | Int  | 15 | Discard any sample that is still running in clockwork rm_contam after this many minutes (set to 0 to never timeout) |  
| timeout_variant_caller | Int  | 120 | Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to 0 to never timeout) |  
| typical_tb_masked_regions | File |  | Bed file of regions to mask when making diff files |  
  
  
## Task-level inputs  
  
### Software settings  
If you are on a backend that does not support call cacheing, you can use the 'bluepeter' inputs to skip the download of the reference genome.  
  
| task | name | type | default | description |  
|:---:|:---:|:---:|:---:|:---:|  
| ClockworkRefPrepTB | bluepeter__tar_indexd_H37Rv_ref | File? |  |  |  
| ClockworkRefPrepTB | bluepeter__tar_indexd_dcontm_ref | File? |  |  |  
| ClockworkRefPrepTB | bluepeter__tar_tb_ref_raw | File? |  |  |  
| cat_reports | out | String  | \'pull_reports.txt\' | Override default output file name with this string |  
| make_mask_and_diff | histograms | Boolean  | false | Should coverage histograms be output? |  
| per_sample_decontam | contam_out_1 | String? |  | Override default output file name with this string |  
| per_sample_decontam | contam_out_2 | String? |  | Override default output file name with this string |  
| per_sample_decontam | counts_out | String? |  | Override default output file name with this string |  
| per_sample_decontam | crash_on_timeout | Boolean  | false | If this task times out, should it stop the whole pipeline (true), or should we just discard this sample and move on (false)? |  
| per_sample_decontam | done_file | String? |  | Override default output file name with this string |  
| per_sample_decontam | no_match_out_1 | String? |  | Override default output file name with this string |  
| per_sample_decontam | no_match_out_2 | String? |  | Override default output file name with this string |  
| per_sample_decontam | subsample_cutoff | Int  | -1 | If a FASTQ file is larger than than size in MB, subsample it with seqtk (set to -1 to disable) |  
| per_sample_decontam | subsample_seed | Int  | 1965 | Seed used for subsampling with seqtk |  
| per_sample_decontam | threads | Int? |  | Try to use this many threads for decontamination. Note that actual number of threads also relies on your hardware. |  
| per_sample_decontam | verbose | Boolean  | true |  |  
| trees | make_nextstrain_subtrees | Boolean  | true |  |  
| trees | outfile | String? |  | Override default output file name with this string |  
| varcall_with_array | crash_on_error | Boolean  | false | If this task, should it stop the whole pipeline (true), or should we just discard this sample and move on (false)? Note that errors that crash the VM (such as running out of space on a GCP instance) will stop the whole pipeline regardless of this setting. |  
| varcall_with_array | crash_on_timeout | Boolean  | false | If this task times out, should it stop the whole pipeline (true), or should we just discard this sample and move on (false)? |  
| varcall_with_array | debug | Boolean  | false | Do not clean up any files and be verbose |  
| varcall_with_array | mem_height | Int? |  | cortex mem_height option. Must match what was used when reference_prepare was run (in other words do not set this variable unless you are also adjusting the reference preparation task) |  
  
  
### Runtime attributes  
These variables adjust runtime attributes, which includes hardware settings. See https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/ for more information.  
  
| task | name | type | default | description |  
|:---:|:---:|:---:|:---:|:---:|  
| cat_reports | disk_size | Int  | 10 | Disk size, in GB. Note that since cannot auto-scale as it cannot anticipate the size of reads from SRA. |  
| get_sample_IDs | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| make_mask_and_diff | addldisk | Int  | 10 | Additional disk size, in GB, on top of auto-scaling disk size. |  
| make_mask_and_diff | cpu | Int  | 8 | Number of CPUs (cores) to request from GCP. |  
| make_mask_and_diff | memory | Int  | 16 | Amount of memory, in GB, to request from GCP. |  
| make_mask_and_diff | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| make_mask_and_diff | retries | Int  | 1 | How many times should we retry this task if it fails after it exhausts all uses of preemptibles? |  
| per_sample_decontam | addldisk | Int  | 100 | Additional disk size, in GB, on top of auto-scaling disk size. |  
| per_sample_decontam | cpu | Int  | 8 | Number of CPUs (cores) to request from GCP. |  
| per_sample_decontam | memory | Int  | 16 | Amount of memory, in GB, to request from GCP. |  
| per_sample_decontam | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| per_sample_decontam | ssd | Boolean  | true | If true, use SSDs for this task instead of HDDs |  
| pull | disk_size | Int  | 100 | Disk size, in GB. Note that since cannot auto-scale as it cannot anticipate the size of reads from SRA. |  
| pull | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| varcall_with_array | addldisk | Int  | 100 | Additional disk size, in GB, on top of auto-scaling disk size. |  
| varcall_with_array | cpu | Int  | 16 | Number of CPUs (cores) to request from GCP. |  
| varcall_with_array | memory | Int  | 32 | Amount of memory, in GB, to request from GCP. |  
| varcall_with_array | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| varcall_with_array | retries | Int  | 1 | How many times should we retry this task if it fails after it exhausts all uses of preemptibles? |  
| varcall_with_array | ssd | Boolean  | true | If true, use SSDs for this task instead of HDDs |  
  

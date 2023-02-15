## Table of Contents  
    * [Workflow-level inputs](#workflow-level-inputs)
    * [Task-level inputs](#task-level-inputs)
      * [Software settings](#software-settings)
      * [Hardware settings](#hardware-settings)
  
  
## Workflow-level inputs  
  
| name | type | default | description |  
|:---:|:---:|:---:|:---:|  
| bad_data_threshold | Float  | 0.05 | If a diff file has higher than this percent (0.5 = 50%) bad data, don not include it in the tree |  
| biosample_accessions | File |  | fastq input -- please see running_myco.md for more information |  
| decorate_tree | Boolean  | false | Should usher, taxonium, and NextStrain trees be generated? Requires input_tree and ref_genome |  
| fastqc_on_timeout | Boolean  | false | If true, fastqc one read from a sample when decontamination times out (see timeout_decontam) |  
| input_tree | File? |  | Base tree to use if decorate_tree = true |  
| min_coverage | Int  | 10 | Positions with coverage below this value will be masked in diff files |  
| paired_fastq_sets | Array |  | fastq input -- please see running_myco.md for more information |  
| ref_genome_for_tree_building | File? |  | Ref genome for building trees -- must have ONLY `&gt;NC_000962.3` on its first line |  
| subsample_cutoff | Int  | 450 | If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable) |  
| subsample_seed | Int  | 1965 | Seed used for subsampling with seqtk |  
| timeout_decontam_part1 | Int  | 20 | Discard any sample that is still running in clockwork map_reads after this many minutes (set to -1 to never timeout) |  
| timeout_decontam_part2 | Int  | 15 | Discard any sample that is still running in clockwork rm_contam after this many minutes (set to -1 to never timeout) |  
| timeout_variant_caller | Int  | 120 | Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to -1 to never timeout) |  
| typical_tb_masked_regions | File |  | Bed file of regions to mask when making diff files |  
  
  
## Task-level inputs  
  
### Software settings  
If you are on a backend that does not support call cacheing, you can use the &#x27;bluepeter&#x27; inputs to skip the download of the reference genome.  
  
| task | name | type | default |  
|:---:|:---:|:---:|:---:|  
| ClockworkRefPrepTB | bluepeter__tar_indexd_H37Rv_ref | File? |  |  
| ClockworkRefPrepTB | bluepeter__tar_indexd_dcontm_ref | File? |  |  
| ClockworkRefPrepTB | bluepeter__tar_tb_ref_raw | File? |  |  
| cat_reports | out | String  | \&quot;pull_reports.txt\&quot; |  
| get_sample_IDs | filter_na | Boolean  | true |  
| make_mask_and_diff | histograms | Boolean  | false |  
| make_mask_and_diff | retries | Int  | 1 |  
| per_sample_decontam | contam_out_1 | String? |  |  
| per_sample_decontam | contam_out_2 | String? |  |  
| per_sample_decontam | counts_out | String? |  |  
| per_sample_decontam | done_file | String? |  |  
| per_sample_decontam | fail_on_timeout | Boolean  | false |  
| per_sample_decontam | filename_metadata_tsv | String  | \&quot;remove_contam_metadata.tsv\&quot; |  
| per_sample_decontam | no_match_out_1 | String? |  |  
| per_sample_decontam | no_match_out_2 | String? |  |  
| per_sample_decontam | subsample_cutoff | Int  | -1 |  
| per_sample_decontam | subsample_seed | Int  | 1965 |  
| per_sample_decontam | threads | Int? |  |  
| per_sample_decontam | verbose | Boolean  | true |  
| pull | fail_on_invalid | Boolean  | false |  
| pull | tar_outputs | Boolean  | false |  
| trees | outfile | String  | \&quot;tree\&quot; |  
| varcall_with_array | debug | Boolean  | false |  
| varcall_with_array | fail_on_timeout | Boolean  | false |  
| varcall_with_array | force | Boolean  | false |  
| varcall_with_array | mem_height | Int? |  |  
| varcall_with_array | retries | Int  | 1 |  
  
  
### Hardware settings  
A note on disk size: On GCP backends, disk size is treated as a maximum. If your task goes above that limit, it will fail.  
  
| task | name | type | default | description |  
|:---:|:---:|:---:|:---:|:---:|  
| cat_reports | disk_size | Int  | 10 | Disk size, in GB. Note that since cannot auto-scale as it cannot anticipate the size of reads from SRA. |  
| get_sample_IDs | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| make_mask_and_diff | addldisk | Int  | 250 | Additional disk size, in GB, on top of auto-scaling disk size. |  
| make_mask_and_diff | cpu | Int  | 16 | Number of CPUs (cores) to request from GCP. |  
| make_mask_and_diff | memory | Int  | 32 | Amount of memory, in GB, to request from GCP. |  
| make_mask_and_diff | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| per_sample_decontam | addldisk | Int  | 100 | Additional disk size, in GB, on top of auto-scaling disk size. |  
| per_sample_decontam | cpu | Int  | 8 | Number of CPUs (cores) to request from GCP. |  
| per_sample_decontam | memory | Int  | 16 | Amount of memory, in GB, to request from GCP. |  
| per_sample_decontam | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| pull | disk_size | Int  | 100 | Disk size, in GB. Note that since cannot auto-scale as it cannot anticipate the size of reads from SRA. |  
| pull | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| varcall_with_array | addldisk | Int  | 250 | Additional disk size, in GB, on top of auto-scaling disk size. |  
| varcall_with_array | cpu | Int  | 16 | Number of CPUs (cores) to request from GCP. |  
| varcall_with_array | memory | Int  | 32 | Amount of memory, in GB, to request from GCP. |  
| varcall_with_array | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
  

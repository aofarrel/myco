## Table of Contents  
    * [Workflow-level inputs](#workflow-level-inputs)
    * [Task-level inputs](#task-level-inputs)
      * [Software settings](#software-settings)
      * [Hardware settings](#hardware-settings)
  
  
## Workflow-level inputs  
  
| name | type | default |  
|:---:|:---:|:---:|  
| bad_data_threshold | Float  | 0.05 |  
| biosample_accessions | File |  |  
| decorate_tree | Boolean  | false |  
| fastqc_on_timeout | Boolean  | false |  
| input_tree | File? |  |  
| min_coverage | Int  | 10 |  
| paired_fastq_sets | Array |  |  
| ref_genome_for_tree_building | File? |  |  
| subsample_cutoff | Int  | 450 |  
| subsample_seed | Int  | 1965 |  
| timeout_decontam_part1 | Int  | 20 |  
| timeout_decontam_part2 | Int  | 15 |  
| timeout_variant_caller | Int  | 120 |  
| typical_tb_masked_regions | File |  |  
  
  
## Task-level inputs  
  
### Software settings  
  
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
  

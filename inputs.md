## Table of Contents  
    * [Workflow-level inputs](#workflow-level-inputs)
    * [Task-level inputs](#task-level-inputs)
      * [Software settings](#software-settings)
      * [Hardware settings](#hardware-settings)
  
  
## Workflow-level inputs  
  
| name | type | default |  
|:---:|:---:|:---:|  
| bad_data_threshold | Float  | default = 0 |  
| biosample_accessions | File |  |  
| decorate_tree | Boolean  | default = false |  
| fastqc_on_timeout | Boolean  | default = false |  
| input_tree | File? |  |  
| min_coverage | Int  | default = 10 |  
| paired_fastq_sets | Array |  |  
| ref_genome_for_tree_building | File? |  |  
| subsample_cutoff | Int  | default = 450 |  
| subsample_seed | Int  | default = 1965 |  
| timeout_decontam_part1 | Int  | default = 20 |  
| timeout_decontam_part2 | Int  | default = 15 |  
| timeout_variant_caller | Int  | default = 120 |  
| typical_tb_masked_regions | File |  |  
  
  
## Task-level inputs  
  
### Software settings  
  
| task | name | type | default |  
|:---:|:---:|:---:|:---:|  
| ClockworkRefPrepTB | bluepeter__tar_indexd_H37Rv_ref | File? |  |  
| ClockworkRefPrepTB | bluepeter__tar_indexd_dcontm_ref | File? |  |  
| ClockworkRefPrepTB | bluepeter__tar_tb_ref_raw | File? |  |  
| cat_reports | out | String  |  |  
| get_sample_IDs | filter_na | Boolean  | default = true |  
| make_mask_and_diff | histograms | Boolean  | default = false |  
| make_mask_and_diff | retries | Int  | default = 1 |  
| per_sample_decontam | contam_out_1 | String? |  |  
| per_sample_decontam | contam_out_2 | String? |  |  
| per_sample_decontam | counts_out | String? |  |  
| per_sample_decontam | done_file | String? |  |  
| per_sample_decontam | fail_on_timeout | Boolean  | default = false |  
| per_sample_decontam | filename_metadata_tsv | String  |  |  
| per_sample_decontam | no_match_out_1 | String? |  |  
| per_sample_decontam | no_match_out_2 | String? |  |  
| per_sample_decontam | subsample_cutoff | Int  |  |  
| per_sample_decontam | subsample_seed | Int  | default = 1965 |  
| per_sample_decontam | threads | Int? |  |  
| per_sample_decontam | verbose | Boolean  | default = true |  
| pull | fail_on_invalid | Boolean  | default = false |  
| pull | tar_outputs | Boolean  | default = false |  
| trees | outfile | String  |  |  
| varcall_with_array | debug | Boolean  | default = false |  
| varcall_with_array | fail_on_timeout | Boolean  | default = false |  
| varcall_with_array | force | Boolean  | default = false |  
| varcall_with_array | mem_height | Int? |  |  
| varcall_with_array | retries | Int  | default = 1 |  
  
  
### Hardware settings  
A note on disk size: On GCP backends, disk size is treated as a maximum. If your task goes above that limit, it will fail.  
  
| task | name | type | default | description |  
|:---:|:---:|:---:|:---:|:---:|  
| cat_reports | disk_size | Int  | default = 10 | Disk size, in GB. Note that since cannot auto-scale as it cannot anticipate the size of reads from SRA. |  
| get_sample_IDs | preempt | Int  | default = 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| make_mask_and_diff | addldisk | Int  | default = 250 | Additional disk size, in GB, on top of auto-scaling disk size. |  
| make_mask_and_diff | cpu | Int  | default = 16 | Number of CPUs (cores) to request from GCP. |  
| make_mask_and_diff | memory | Int  | default = 32 | Amount of memory, in GB, to request from GCP. |  
| make_mask_and_diff | preempt | Int  | default = 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| per_sample_decontam | addldisk | Int  | default = 100 | Additional disk size, in GB, on top of auto-scaling disk size. |  
| per_sample_decontam | cpu | Int  | default = 8 | Number of CPUs (cores) to request from GCP. |  
| per_sample_decontam | memory | Int  | default = 16 | Amount of memory, in GB, to request from GCP. |  
| per_sample_decontam | preempt | Int  | default = 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| pull | disk_size | Int  | default = 100 | Disk size, in GB. Note that since cannot auto-scale as it cannot anticipate the size of reads from SRA. |  
| pull | preempt | Int  | default = 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| varcall_with_array | addldisk | Int  | default = 250 | Additional disk size, in GB, on top of auto-scaling disk size. |  
| varcall_with_array | cpu | Int  | default = 16 | Number of CPUs (cores) to request from GCP. |  
| varcall_with_array | memory | Int  | default = 32 | Amount of memory, in GB, to request from GCP. |  
| varcall_with_array | preempt | Int  | default = 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
  

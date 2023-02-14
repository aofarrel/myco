## Table of Contents  
    * [Workflow-level inputs](#workflow-level-inputs)
    * [Task-level inputs](#task-level-inputs)
      * [Software settings](#software-settings)
      * [Hardware settings](#hardware-settings)
  
  
## Workflow-level inputs  
  
| name | type | default |  
|:---:|:---:|:---:|  
| biosample_accessions | File |  |  
| less_scattering | Boolean  | default = false |  
| min_coverage | Int  | default = 10 |  
| paired_fastqs | Array |  |  
| subsample_cutoff | Int  |  |  
| subsample_seed | Int  | default = 1965 |  
| tar_fqs | Boolean  | default = false |  
| typical_tb_masked_regions | File |  |  
  
  
## Task-level inputs  
  
### Software settings  
  
| task | name | type | default |  
|:---:|:---:|:---:|:---:|  
| ClockworkRefPrepTB | bluepeter__tar_indexd_H37Rv_ref | File? |  |  
| ClockworkRefPrepTB | bluepeter__tar_indexd_dcontm_ref | File? |  |  
| ClockworkRefPrepTB | bluepeter__tar_tb_ref_raw | File? |  |  
| cat_reports | keep_only_unique_lines | Boolean? | default = false |  
| cat_reports | out_filename | String? |  |  
| decontaminate_many_samples | contam_out_1 | String? |  |  
| decontaminate_many_samples | contam_out_2 | String? |  |  
| decontaminate_many_samples | counts_out | String? |  |  
| decontaminate_many_samples | done_file | String? |  |  
| decontaminate_many_samples | filename_metadata_tsv | String  |  |  
| decontaminate_many_samples | no_match_out_1 | String? |  |  
| decontaminate_many_samples | no_match_out_2 | String? |  |  
| decontaminate_many_samples | threads | Int? |  |  
| decontaminate_many_samples | verbose | Boolean  | default = true |  
| decontaminate_one_sample | contam_out_1 | String? |  |  
| decontaminate_one_sample | contam_out_2 | String? |  |  
| decontaminate_one_sample | counts_out | String? |  |  
| decontaminate_one_sample | done_file | String? |  |  
| decontaminate_one_sample | filename_metadata_tsv | String  |  |  
| decontaminate_one_sample | no_match_out_1 | String? |  |  
| decontaminate_one_sample | no_match_out_2 | String? |  |  
| decontaminate_one_sample | threads | Int? |  |  
| decontaminate_one_sample | verbose | Boolean  | default = true |  
| make_mask_and_diff | retries | Int  | default = 1 |  
| make_mask_and_diff_ | retries | Int  | default = 1 |  
| pull | disk_size | Int  | default = 60 |  
| pull | fail_on_invalid | Boolean  | default = false |  
| pull | subsample_cutoff | Int  | default = 450 |  
| pull | subsample_seed | Int  | default = 1965 |  
| varcall_with_array | debug | Boolean  | default = false |  
| varcall_with_array | force | Boolean  | default = false |  
| varcall_with_array | mem_height | Int? |  |  
| varcall_with_array | retries | Int  | default = 1 |  
| varcall_with_tarballs | debug | Boolean  | default = true |  
| varcall_with_tarballs | force | Boolean  | default = false |  
| varcall_with_tarballs | mem_height | Int? |  |  
| varcall_with_tarballs | reads_files | Array |  |  
| varcall_with_tarballs | retries | Int  | default = 1 |  
| varcall_with_tarballs | sample_name | String? |  |  
  
  
### Hardware settings  
  
| task | name | type | default |  
|:---:|:---:|:---:|:---:|  
| cat_reports | preempt | Int? | default = 1 |  
| decontaminate_many_samples | addldisk | Int  | default = 100 |  
| decontaminate_many_samples | cpu | Int  | default = 16 |  
| decontaminate_many_samples | memory | Int  | default = 32 |  
| decontaminate_many_samples | preempt | Int  | default = 0 |  
| decontaminate_one_sample | addldisk | Int  | default = 100 |  
| decontaminate_one_sample | cpu | Int  | default = 8 |  
| decontaminate_one_sample | memory | Int  | default = 16 |  
| decontaminate_one_sample | preempt | Int  | default = 1 |  
| get_sample_IDs | preempt | Int? | default = 1 |  
| make_mask_and_diff | addldisk | Int  | default = 250 |  
| make_mask_and_diff | cpu | Int  | default = 16 |  
| make_mask_and_diff | memory | Int  | default = 32 |  
| make_mask_and_diff | preempt | Int  | default = 1 |  
| make_mask_and_diff_ | addldisk | Int  | default = 250 |  
| make_mask_and_diff_ | cpu | Int  | default = 16 |  
| make_mask_and_diff_ | memory | Int  | default = 32 |  
| make_mask_and_diff_ | preempt | Int  | default = 1 |  
| pull | preempt | Int  | default = 1 |  
| varcall_with_array | addldisk | Int  | default = 250 |  
| varcall_with_array | cpu | Int  | default = 16 |  
| varcall_with_array | memory | Int  | default = 32 |  
| varcall_with_array | preempt | Int  | default = 1 |  
| varcall_with_tarballs | addldisk | Int  | default = 250 |  
| varcall_with_tarballs | cpu | Int  | default = 16 |  
| varcall_with_tarballs | memory | Int  | default = 32 |  
| varcall_with_tarballs | preempt | Int  | default = 1 |  
  

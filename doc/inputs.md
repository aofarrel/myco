## Table of Contents  
    * [Workflow-level inputs](#workflow-level-inputs)
    * [Task-level inputs](#task-level-inputs)
      * [Software settings](#software-settings)
      * [Runtime attributes](#runtime-attributes)
  
See /inputs/example_inputs.json for examples.  
  
## Workflow-level inputs  
  
| name | type | default | description |  
|:---:|:---:|:---:|:---:|  
| bad_data_threshold | Float  | 0.05 | If a diff file has higher than this percent (0.5 = 50%) bad data, do not include it in the tree |  
| biosample_accessions | File |  | fastq input -- please see running_myco.md for more information |  
| decorate_tree | Boolean  | false | Should usher, taxonium, and NextStrain trees be generated? Requires input_tree and ref_genome |  
| fastqc_on_timeout | Boolean  | false | If true, fastqc one read from a sample when decontamination or variant calling times out |  
| force_diff | Boolean  | false | If true and if decorate_tree is false, generate diff files. (Diff files will always be created if decorate_tree is true.) |  
| input_tree | File? |  | Base tree to use if decorate_tree = true |  
| min_coverage | Int  | 10 | Positions with coverage below this value will be masked in diff files |  
| paired_fastq_sets | Array |  | fastq input -- please see running_myco.md for more information |  
| ref_genome_for_tree_building | File? |  | Ref genome for building trees -- must have ONLY `>NC_000962.3` on its first line |  
| subsample_cutoff | Int  | 450 | If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable) |  
| subsample_seed | Int  | 1965 | Seed used for subsampling with seqtk |  
| timeout_decontam_part1 | Int  | 20 | Discard any sample that is still running in clockwork map_reads after this many minutes (set to -1 to never timeout) |  
| timeout_decontam_part2 | Int  | 15 | Discard any sample that is still running in clockwork rm_contam after this many minutes (set to -1 to never timeout) |  
| timeout_variant_caller | Int  | 120 | Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to -1 to never timeout) |  
| typical_tb_masked_regions | File |  | Bed file of regions to mask when making diff files |  
  
  
## Task-level inputs  
  
### Software settings  
If you are on a backend that does not support call cacheing, you can use the `bluepeter` inputs to skip the download of the reference genome.  
  
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
| per_sample_decontam | subsample_cutoff | Int  | -1 | If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable) |  
| per_sample_decontam | subsample_seed | Int  | 1965 | Seed used for subsampling with seqtk |  
| per_sample_decontam | threads | Int? |  | Try to use this many threads for decontamination. Note that actual number of threads also relies on your hardware. |  
| per_sample_decontam | verbose | Boolean  | true |  |  
| trees | outfile | String  | \'tree\' | Override default output file name with this string |  
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
  

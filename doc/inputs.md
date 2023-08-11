### A quick introduction to input variables in WDL
There are two kinds of input user-accessible variables in WDL: Workflow-level inputs and task-level inputs. If you are using Terra, you probably don't need to know anything about the difference between them except that task-level inputs get alphabetically sorted in Terra's UI below workflow-level inputs. 

This pipeline uses a lot of external tools, and I tend to WDLize every possible input variable, so there are a lot of input variables in this pipeline. **The vast majority of them are optional.** What's most important is your fastqs.

# Workflow-level inputs  
  
## FASTQs  
Each version of myco has a slightly different way of inputting FASTQs. A basic explanation for each workflow is in the table below. You can find more detailed explanations in each workflow's workflow-level readme.  
  
| name | type | workflow | description |  
|:---:|:---:|:---:|:---:|  
| biosample_accessions | File | myco_sra | File of BioSample accessions to pull, one accession per line |  
| paired_decontaminated_fastq_sets | Array | myco_cleaned | Nested array of decontaminated and merged fastq pairs. Each inner array represents one sample; each sample needs precisely one forward read and one reverse read. |  
| paired_fastq_sets | Array | myco_raw | Nested array of paired fastqs, each inner array representing one samples worth of paired fastqs |  
  
Regardless of which version of myco you use, please make sure your FASTQs:
* is Illumina paired-end data <sup>†</sup>  
* is grouped per-sample   
* len(quality scores) = len(nucleotides) for every line <sup>†</sup>  
* is actually [MTBC](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=77643)  
<sup>†</sup> myco_sra.wdl is able to detect these issues and will throw out those samples without erroring. Other forms of myco are not able to detect these issues.
It is recommend that you also keep an eye on the total size of your FASTQs. Individual files over subsample_cutoff (default450 MB, -1 disables this check) will be downsampled, but keep an eye on the cumulative size of samples. For example, a sample like SAMEA968096 has 12 run accessions associated with it. Individually, none of these run accessions' FASTQs are over 1 GB in size, but the sum total of these FASTQs could quickly fill up your disk space. (You probably should not be using SAMEA968096 anyway because it is in sample group, which can cause other issues.)

myco_cleaned expects that the FASTQs you are putting into have already been cleaned and merged. It's recommend you do this by running [Decontam_and_Combine](https://dockstore.org/workflows/github.com/aofarrel/clockwork-wdl/Decontam_And_Combine_One_Samples_Fastqs).  
 
## Quality control
| name | type | myco_sra default | description |  
|:---:|:---:|:---:|:---:|  
| covstats_qc_cutoff_coverages | Float  | 10 | If covstats thinks coverage is below this, throw out this sample |  
| covstats_qc_cutoff_unmapped | Float  | 2 | If covstats thinks this percentage (50 = 50%) of data does not map to H37Rv, throw out this sample |  
| covstats_qc_skip_entirely | Boolean  | false | Should we avoid running covstats? Does not affect other forms of QC. |  
| diff_max_pct_low_coverage | Float  | 0.05 (myco_sra), 0.20 (myco_raw) | If more than this percent (0.5 = 50%) of a sample's sites get masked for being below `diff_min_coverage_per_site`, throw out the whole sample. |  
| diff_min_coverage_per_site | Int  | 10 | Positions with coverage below this value will be masked in diff files |
| early_qc_apply_cutoffs | Boolean  | false | If true, run fastp + TBProfiler on decontaminated fastqs and apply cutoffs to determine which samples should be thrown out. |  
| early_qc_cutoff_q30 | Float  | 0.9 | Decontaminated samples with less than this percentage (as float, 0.5 = 50%) of reads above qual score of 30 will be discarded iff early_qc_apply_cutoffs is also true. |  
| early_qc_skip_entirely | Boolean  | false | Do not run early QC (fastp + fastq-TBProfiler) at all. Does not affect whether or not TBProfiler is later run on bams. Overrides early_qc_apply_cutoffs. |  
| fastqc_on_timeout | Boolean  | false | (myco_sra only) If true, fastqc one read from a sample when decontamination or variant calling times out |  

Note that all forms of QC will throw out entire samples, with two exceptions: 
  * `diff_min_coverage_per_site` works on a per-site basis, masking individual sites below that value but not throwing out the entire sample (unless so many sites get masked that `diff_max_pct_low_coverage` kicks in)
  * `fastqc_on_timeout` will run fastQC if true, but does not parse its output - the samples have already been discarded upstream when they timed out
  
## Timeouts  
When working with data of unknown quality, it can be helpful to quickly remove samples that are likely low-quality. While developing myco on SRA data, we noticed that if a given sample took an unusually long time in the decontamination or variant calling step, they were likely to end up filtered out by the final quality control steps of the pipeline. This is especially true of the decontamination step -- the more contamination a sample has, the more that step has to do. This heuristic was defined on the default runtime attributes and using Terra as a backend, so straying from those defaults is likely to make the default timeout values less useful. This *includes* changing from SDDs to HDDs! 
  
| name | type | myco_sra default | description |  
|:---:|:---:|:---:|:---:|  
| timeout_decontam_part1 | Int  | 20<sup>†</sup> | Discard any sample that is still running in clockwork map_reads after this many minutes (set to 0 to never timeout |
| timeout_decontam_part2 | Int  | 15<sup>†</sup>  | Discard any sample that is still running in clockwork rm_contam after this many minutes (set to 0 to never timeout) |  
| timeout_variant_caller | Int  | 120<sup>†</sup> | Discard any sample that is still running in clockwork variant_call_one_sample after this many minutes (set to 0 to never timeout) | 

<sup>†</sup> myco_raw and myco_cleaned default to not using this heuristic at all, so their defaults are 0.


## Variant caller inputs
Usually, I write WDLs in a way that makes their runtime attributes and rarely-used optional arguments task-level, and everything else workflow-level. However, myco uses some workarounds that require it to technically have three copies of the variant caller task, which means that if I didn't make the variant caller's inputs workflow-level, there would be three sets of task-level inputs for the variant caller.

| name | type | myco_sra default | description |  
|:---:|:---:|:---:|:---:|  
| variantcalling_addl_disk        | Int     | 100   | Additional disk size, in GB, on top of auto-scaling disk size. |
| variantcalling_cpu              | Int     | 16    | Number of CPUs (cores) to request from GCP. |
| variantcalling_crash_on_error   | Boolean | false | If this task errors out, should it stop the whole pipeline (true), or should we just discard this sample and move on (false)? Note that errors that crash the VM (such as running out of space on a GCP instance) will stop the whole pipeline regardless of this setting. |  
| variantcalling_crash_on_timeout | Boolean | false | If this task times out, should it stop the whole pipeline (true), or should we just discard this sample and move on (false)? | 
| variantcalling_debug            | Boolean | false | Do not clean up any files and be verbose |  
| variantcalling_mem_height       | Int?    |       | cortex mem_height option. Must match what was used when reference_prepare was run (in other words do not set this variable unless you are also adjusting the reference preparation task) | 
| variantcalling_memory           | Int     | 32    | Amount of memory, in GB, to request from GCP. |
| variantcalling_preemptibles     | Int     |  1    | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| variantcalling_retries          | Int     |  1    | How many times should we retry this task if it fails after it exhausts all uses of preemptibles? | 
| variantcalling_ssd              | Boolean | true  | If true, use SSDs for this task instead of HDDs  | 

  
## Miscellanous workflow-level inputs  
  
| name | type | default | description |  
|:---:|:---:|:---:|:---:|  
| diff_force | Boolean  | false | If true and if decorate_tree is false, generate diff files. (Diff files will always be created if decorate_tree is true.) |  
| diff_mask_these_regions | File? | [this CRyPTIC mask file](https://github.com/iqbal-lab-org/cryptic_tb_callable_mask/blob/44f884558bea4ee092ce7c5c878561200fcee92f/R00000039_repregions.bed) | Bed file of regions to mask when making diff files |  
| quick_tasks_disk_size | Int  | 10 | If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)isk size in GB to use for quick file-processing tasks; increasing this might slightly speed up file localization |  
| subsample_cutoff | Int  | 450 | If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable) |  
| subsample_seed | Int  | 1965 | Seed used for subsampling with seqtk |  
| tbprofiler_on_bam | Boolean  | varies<sup>†</sup> | If true, run TBProfiler on BAMs |  
| tree_decoration | Boolean  | false | Should usher, taxonium, and NextStrain trees be generated? |  
| tree_to_decorate | File? | [this draft tree](https://console.cloud.google.com/storage/browser/_details/topmed_workflow_testing/tb/trees/alldiffs_mask2ref.L.fixed.pb;tab=live_object) | Base tree to use if decorate_tree = true |  

<sup>†</sup> Defaults to true for myco_sra and false for myco_raw for historical reasons
  

  
# Task-level inputs 

Many of these settings just change the name of output files or runtime attributes, but for the sake of ease-of-use we wanted to include them on a single list in the alphabetical order they show up in on Terra. As such, settings that DO NOT relate to filename outputs nor runtime attributes are in bold.

For more info on runtime settings, see https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/


| task | name | type | default | description |  
|:---:|:---:|:---:|:---:|:---:|  
| **check_fastqs** | **output_fastps_cleaned_fastqs** | **Boolean**  | **false** | **Use fastp's cleaned fastqs for subsequent tasks (true) instead of discarding them (false). Setting this to false means you are effectively only using fastp to check if a sample is valid, keeping or throwing out the entire sample based on this information.** |  
| collate_depth | disk_size | Int  | 10 | Disk size, in GB. This task cannot autoscale as it cannot anticipate the size of reads from SRA.  |  
| collate_resistance | disk_size | Int  | 10 | Disk size, in GB. This task cannot autoscale as it cannot anticipate the size of reads from SRA.  |  
| collate_strains | disk_size | Int  | 10 | Disk size, in GB. This task cannot autoscale as it cannot anticipate the size of reads from SRA. |  
| decontam_each_sample | addldisk | Int  | 100 | Additional disk size, in GB, on top of auto-scaling disk size. |  
| decontam_each_sample | cpu | Int  | 8 | Number of CPUs (cores) to request from GCP. |  
| decontam_each_sample | contam_out_1 | String? |  |  |  
| decontam_each_sample | contam_out_2 | String? |  |  |  
| decontam_each_sample | counts_out | String? |  |  |  
| **decontam_each_sample** | **crash_on_timeout** | **Boolean**  | **false** | **If this task times out, should it stop the whole pipeline (true), or should we just discard this sample and move on (false)?** |  
| decontam_each_sample | done_file | String? |  |  |  
| decontam_each_sample | memory | Int  | 16 | Amount of memory, in GB, to request from GCP. |  
| decontam_each_sample | no_match_out_1 | String? |  | |  
| decontam_each_sample | no_match_out_2 | String? |  |  |  
| decontam_each_sample | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| decontam_each_sample | ssd | Boolean  | true | If true, use SSDs for this task instead of HDDs |  
| **decontam_each_sample** | **subsample_cutoff** | **Int**  | **-1** | **If a FASTQ file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)** |  
| **decontam_each_sample** | **subsample_seed**   | **Int**  | **1965** | **Seed used for subsampling with seqtk** |  
| **decontam_each_sample** | **threads**          | **Int?** |  | **Try to use this many threads for decontamination. Note that actual number of threads also relies on your hardware.** |  
| **decontam_each_sample** | **verbose**          | **Boolean**  | **true** | **Enable/Disable debug mode** |  
| **FastqcWF** | **limits** | **File?** |  | **Limits file defining fastQC cutoffs** |  
| get_sample_IDs | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| make_mask_and_diff | addldisk | Int  | 10 | Additional disk size, in GB, on top of auto-scaling disk size. |  
| make_mask_and_diff | cpu | Int  | 8 | Number of CPUs (cores) to request from GCP. |  
| **make_mask_and_diff** | **histograms** | **Boolean**  | **false** | **Should coverage histograms be output?** |  
| make_mask_and_diff | memory | Int  | 16 | Amount of memory, in GB, to request from GCP. |  
| make_mask_and_diff | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| make_mask_and_diff | retries | Int  | 1 | How many times should we retry this task if it fails after it exhausts all uses of preemptibles? |  
| merge_reports | disk_size | Int  | 10 | Disk size, in GB. This task cannot autoscale as it cannot anticipate the size of reads from SRA. |  
| profile | bam_suffix | String? |  |  |  
| profile | cpu | Int  | 2 | Number of CPUs (cores) to request from GCP. |  
| profile | memory | Int  | 4 | Amount of memory, in GB, to request from GCP. |  
| profile | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| profile | ssd | Boolean  | false | If true, use SSDs for this task instead of HDDs |  
| pull | disk_size | Int  | 100 | Disk size, in GB. This task cannot autoscale as it cannot anticipate the size of reads from SRA. |  
| pull | preempt | Int  | 1 | How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance? |  
| **trees** | **detailed_clades** | **Boolean**  | **false** | **If true, run usher sampled diff with `-D` flag** |  
| trees | in_prefix_summary | String? |  |  |  
| **trees** | **make_nextstrain_subtrees** | **Boolean**  | **true** | **Should our Nextstrain (Auspice-compatiable) output tree consist of multiple subtrees (true), or just one big tree (false)? Setting this to true might be useful if you have trouble loading very large trees into Auspice.** |  
| **trees** | **metadata_tsv** | **File?** |  | **Metadata TSV to annotate the input tree with** |  
| trees | out_diffs | String  | \'_combined\' |  |  
| trees | out_prefix | String  | \'tree\' |  |  
| trees | out_prefix_summary | String? |  |  |  
| trees | out_tree_annotated_pb | String  | \'_annotated\' |  |  
| trees | out_tree_nextstrain | String  | \'_auspice\' |  |  
| trees | out_tree_nwk | String  | \'_nwk\' |  |  
| trees | out_tree_raw_pb | String  | \'_raw\' |  |  
| trees | out_tree_taxonium | String  | \'_taxonium\' |  |  
| **trees** | **reroot_to_this_node** | **String?** |  | **ID of the node you want to reroot to. Usually, this ID is a BioSample accession.** |  
| **trees** | **subtree_only_new_samples** | **Boolean**  | **true** | **Iff make_nextstrain_subtrees is true, the subtrees are focused only on newly-added samples (eg, samples that went through this run of the pipeline, rather than samples that were already on the pre-existing tree)** |  
| **trees** | **summarize_input_mat** | **Boolean**  | **true** | **Run `matutils summary` on the input tree before adding new samples to it, and output that summary** |  

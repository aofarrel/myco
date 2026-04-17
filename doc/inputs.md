### A quick introduction to input variables in WDL
There are two kinds of user-accessible input variables in WDL: Workflow-level inputs and task-level inputs. If you are using Terra, you probably don't need to know anything about the difference between them except that task-level inputs get alphabetically sorted in Terra's UI below workflow-level inputs. 

This pipeline uses a lot of external tools, and I tend to WDLize every possible input variable, so there are a lot of input variables in this pipeline. **The vast majority of them are optional.** What's most important is your fastqs.

# Workflow-level inputs  
  
## FASTQs  
Each version of myco has a slightly different way of inputting FASTQs. A basic explanation for each workflow is in the table below. You can find more detailed explanations in each workflow's workflow-level readme.  
  
| name | type | workflow | description |  
|:---:|:---:|:---:|:---:|  
| biosample_accessions | File | myco_sra | File of BioSample accessions to pull, one accession per line |  
| paired_decontaminated_fastq_sets | Array | myco_simple| Nested array of decontaminated and merged fastq pairs. Each inner array represents one sample; each sample needs precisely one forward read and one reverse read. |  
| paired_fastq_sets | Array | myco_raw | Nested array of paired fastqs, each inner array representing one samples worth of paired fastqs |  
  
Regardless of which version of myco you use, please make sure your FASTQs:
* is Illumina paired-end data <sup>†</sup>  
* is grouped per-sample   
* len(quality scores) = len(nucleotides) for every line <sup>†</sup>  
* is actually [MTBC](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=77643)  

<sup>†</sup> myco_sra.wdl is able to detect these issues and will throw out those samples without erroring. Other forms of myco are not able to detect these issues.  

It is recommend that you also keep an eye on the total size of your FASTQs. Individual files over subsample_cutoff (default450 MB, -1 disables this check) will be downsampled, but keep an eye on the cumulative size of samples. For example, a sample like SAMEA968096 has 12 run accessions associated with it. Individually, none of these run accessions' FASTQs are over 1 GB in size, but the sum total of these FASTQs could quickly fill up your disk space. (You probably should not be using SAMEA968096 anyway because it is in sample group, which can cause other issues.)

myco_simple expects that the FASTQs you are putting into have already been cleaned and decontaminated, but this isn't a hard requirement. What is a hard requirement is that you have precisely one forward read and one reverse read per sample -- if you have multi-lane samples across various fastqs, they will need to be merged first.
 
## Other inputs
*This list is incomplete. Please see parameter_meta for more information.*  

| name | type | myco_raw default | myco_sra default | description |  
|:---:|:---:|:---:|:---:|  
| comment | String? |  |  | Passed directly as a workflow output `tbd_comment`, useful for Terra data tables in some scenarios |  
| call_as_reference_bedfile | File? | [this CRyPTIC mask file](https://github.com/iqbal-lab-org/cryptic_tb_callable_mask/blob/44f884558bea4ee092ce7c5c878561200fcee92f/R00000039_repregions.bed) | same as myco_raw | Bed file of regions to mask as reference when making diff files |  
| guardrail_mode | Boolean | true | true | Implements safeguards, see section below for more information | 
| sample_max_pct_masked | Int | | Samples who have more than this percent (as int, 50 = 50%) of positions with coverage below site_min_depth will be discarded |
| sample_min_q30 | Int | 80 | 90 | Decontaminated samples with less than this percent (as int, 50 = 50%) of reads above qual score of 30 will be discarded |
| skip_covstats | Boolean  | true | true | Skip covstats |  
| subsample_cutoff | Int  | -1 (disabled) | 450 (ie, 450 MB) | If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable) |  
| subsample_reads | Int | 2000000 | 1000000 | If a fastq file is larger than `subsample_cutoff`, downsample it to this many reads |

> [!WARNING]  
> If `just_like_2024` is set to non-default value true, `subsample_cutoff` will be overwritten to 450 (as in 450 megabytes), and `subsample_reads` to 1000000 (as in 1,000,000 reads). These values were chosen carefully to help balance for the high cost of running myco_sra on the entirity of NCBI SRA's tuberculosis data while retaining suitable depth, but may be unacceptably small for other use cases.


> [!NOTE]  
> Regions within call_as_reference_bedfile are called as reference in resulting diff file. In all other scenarios (indels, low coverage, ambigious call) when we refer to a "masked" position, we mean one that is explictly included in the diff file as `-`. There is a minor distinction between "reference" and "masked" with regard to placement of samples on a phylogenetic tree by UShER/[Tree Nine](https://github.com/aofarrel/tree_nine).
>
> Be aware that `call_as_reference_bedfile` only applies to the diff file. The VCF is not affected and may call non-referennce variants in those regions.


## Guardrail Mode  
Guardrail Mode implements timers to certain myco tasks, which help prevent edge case samples from causing runaway cloud costs and pipeline stalling. This is especially important in the decontamination step, as decontamination occurs before most QC checks and requires a long time to complete on extremely contaminated large samples, which are ultimately doomed to fail QC checks anyway. The defaults myco uses for guardrails are relatively lenient to ensure the maximum number of likely-to-pass samples make it through the pipeline. It's recommended to leave this enabled unless your fastqs are huge, or you are running on slow HDDs.

In previous versions, `subsample_cutoff` was in some cases overwritten by a default if `guardrail_mode` = True. In this version, this is no longer the case. subsample_cutoff will always be respected regardless of the value of guardrail_mode. That being said, it is one of the best possible guardrails against bad data, so consider leaving it to the default value.


## Deprecated

decontam_use_CDC_varpipe_ref: "If true, use CDC varpipe decontamination reference. If false, use CRyPTIC decontamination reference."  
CDC uses their own version of clockwork's decontamination reference, which I call "CDC varpipe" since I pulled it from the varpipe repo. This is currently a null op in myco_raw as it'd require I maintain double the number of Docker images, and it doesn't delineate between human vs NTM vs other forms of contamination (which the task currently requires for some outputs). If there is a demand for CDC varpipe I can make this an option again, but I nevertheless gently recommend against using it due to unclear provenance and contents.

subsample_seed: "Seed used for subsampling with seqtk"  
To reduce the number of inputs, this is now 1965.  

tbprofiler_on_bam: "If true, run TBProfiler on BAMs"


  
# Task-level inputs 
Many of myco's tasks have exposed runtime attributes. These don't show up on the workflow page, but do show up in Terra, and can be accessed by any WDL executor by referencing them in the input JSON. Most of them are for configuring runtime attributes: https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/


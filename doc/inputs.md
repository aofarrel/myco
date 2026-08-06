### A quick introduction to input variables in WDL
There are two kinds of user-accessible input variables in WDL: Workflow-level inputs and task-level inputs. If you are using Terra, you probably don't need to know anything about the difference between them except that task-level inputs get alphabetically sorted in Terra's UI below workflow-level inputs. 

I tend to expose a lot of variables when writing pipelines, so there a lot of possible input variables in myco. **The vast majority of them are optional.**

# Workflow-level inputs  
  
## FASTQs  
Each version of myco has a slightly different way of inputting FASTQs. A basic explanation for each workflow is in the table below. You can find more detailed explanations in each workflow's workflow-level readme.  
  
| name | type | workflow | description |  
|:---:|:---:|:---:|:---:|  
| biosample_accessions | File | myco_sra | File of BioSample accessions to pull, one accession per line |  
| paired_decontaminated_fastq_sets | Array | myco_simple| Nested array of decontaminated and merged fastq pairs. Each inner array represents one sample; each sample needs precisely one forward read and one reverse read. |  
| paired_fastq_sets | Array | myco_raw | Nested array of paired fastqs, each inner array representing one samples worth of paired fastqs |  
  
Regardless of which version of myco you use, please make sure your FASTQs:
* are Illumina paired-end data <sup>†</sup>  
* are grouped per-sample   
* len(quality scores) = len(nucleotides) for every line <sup>†</sup>  
* are actually [MTBC](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=77643)  

<sup>†</sup> myco_sra.wdl is able to detect these issues and will throw out those samples without erroring. Other forms of myco are not able to detect these issues.  

myco_simple additionally expects that the FASTQs you are putting into have already been cleaned and decontaminated, but this isn't a hard requirement. What is a hard requirement is that you have precisely one forward read and one reverse read per sample -- if you have multi-lane samples across various fastqs, they will need to be merged first.

It is recommend that you also keep an eye on the total size of your FASTQs. Individual files over subsample_cutoff (default 450 MB, -1 disables this check) will be downsampled, but keep an eye on the cumulative size of samples. For example, if a BioSample has 12 run accessions associated with, it's possible none of those run accessions FASTQs are very large, but the sum total of these FASTQs may be more than you bargained for.

### Here be dragons: Beware of broken samples on NCBI SRA
Some samples on NCBI SRA are corrupted, improperly tagged, "broken", very large, etc. In the name of running on as many samples as possible with minimal curation, myco_sra relies on [SRANWRP](https://github.com/aofarrel/SRANWRP/), which itself relies on an intentionally outdated version of fasterq-dump to pull samples from SRA in order to only allow Illumina samples to be downloaded (as clockwork can only handle Illumina PE samples). This can handle BioSamples of mixed type; if a BioSample has Illumina and PacBio reads, only the Illumnina will be pulled and processed by myco_sra, and this will occur without throwing an error. SRANWRP additionally accounted for almost all other edge cases within the entire collection of MTBC samples that I was aware of as of mid-2023, such as samples with corrupted quality scores.

An exception: A small number of BioSamples such as [SAMEA968096](https://www.ncbi.nlm.nih.gov/sra/?term=SAMEA968096) are in "sample pools." Such samples will return all run accessions for all BioSamples in that sample pool, often including a generic barcode sample. This may cause issues with downstream analysis; such samples should likely be avoided. There does not appear to be a consistent way to detect sample groups using fasterq-dump, which is why this currently is not handled automatically, but samples within a sample pool are marked on SRA's UI: <img width="382" height="26" alt="Screenshot 2026-08-06 at 3 21 58 PM" src="https://github.com/user-attachments/assets/2622ba59-3040-4492-a442-4638f83fee8f" />

[I have a compiled denylists of all problematic MTBC samples that I'm aware of](https://github.com/aofarrel/SRANWRP/tree/main/inputs/denylists), including samples within sample pools, but there may be more out there.
 
## Other inputs
*This is just a brief overview. Please see myco_sra and myco_raw's respective parameter_meta sections within their WDLs for more information.*  

| name | type | myco_raw default | myco_sra default | description |  
|:---:|:---:|:---:|:---:|:---:|    
| call_as_reference_bedfile | File? | [this CRyPTIC mask file](https://github.com/iqbal-lab-org/cryptic_tb_callable_mask/blob/44f884558bea4ee092ce7c5c878561200fcee92f/R00000039_repregions.bed) | same as myco_raw | Bed file of regions to mask as reference when making diff files **(be aware the VCF is not affected by this file, only the final MAPLE diff)** |  
| comment | String? |  |  | Passed directly as a workflow output `tbd_comment`, useful for Terra data tables in some scenarios |  
| guardrail_mode | Boolean | true | true | Implements safeguards, see section below for more information | 
| just_like_2024 | Boolean | false | false | Override a bunch of QC metrics, as well as `subsample_cutoff` and `subsample_reads` [to match what was used to build UCSC's SRA tree](10.1101/2025.07.22.25331806) **(do not enable this unless you know what you're doing; this exists in the name of publication reproducibility)** | 
| sample_max_pct_masked | Int | 20 | 20 | Samples who have more than this percent (as int, so 20 = 20%) of positions with coverage below site_min_depth will be discarded |
| sample_min_q30 | Int | 80 | 80 | Decontaminated samples with less than this percent (as int, so 80 = 80%) of reads above qual score of 30 will be discarded (intentionally lower than CDC's preferred 85% as in our experience 85% cutoff removes too CalTBNet samples) |
| site_min_depth | Int | 10 | 10 | Positions with coverage below this value will be explicitly masked (`-`) in diff files; see also sample_max_pct_masked |  
| subsample_cutoff | Int  | -1 (disabled) | 450 | If a fastq file is larger than than size in MB (megabytes, not megabases), subsample it with seqtk (set to -1 to disable) to `subsample_reads` number of reads |  
| subsample_reads | Int | 2000000 | 1000000 | If a fastq file is larger than `subsample_cutoff`, downsample it to this many reads |

> [!WARNING]  
> If `just_like_2024` is set to non-default value true, in additional to several QC variables being changed, `subsample_cutoff` will be overwritten to 450 (as in 450 megabytes), and `subsample_reads` to 1000000 (as in 1,000,000 reads). These values were chosen carefully to help balance for the high cost of running myco_sra on the entirety of NCBI SRA's tuberculosis data while retaining suitable depth, but may be unacceptably small for other use cases.

> [!NOTE]  
> Regions within call_as_reference_bedfile are called as reference in resulting diff file, which is to say, they are not mentioned at all (in the same way a VCF doesn't mention every site that is reference). In all other scenarios (indels, low coverage as per `site_min_depth`, ambiguous call) when we refer to a "masked" position, we mean one that is explicitly included in the diff file as `-`. There is a minor distinction between "reference" and "masked" with regard to placement of samples on a phylogenetic tree by UShER/[Tree Nine](https://github.com/aofarrel/tree_nine).

## Guardrail Mode  
Guardrail Mode implements timers to certain myco tasks, which help prevent edge case samples from causing runaway cloud costs and pipeline stalling. This is especially important in the decontamination step, as decontamination occurs before most QC checks and requires a long time to complete on extremely contaminated large samples, which are ultimately doomed to fail QC checks anyway. The defaults myco uses for guardrails are relatively lenient to ensure the maximum number of likely-to-pass samples make it through the pipeline. It's recommended to leave this enabled unless your fastqs are huge, or you are running on slow HDDs.

In previous versions, `subsample_cutoff` was in some cases overwritten by a default if `guardrail_mode` = True. In this version, this is no longer the case. subsample_cutoff will always be respected regardless of the value of guardrail_mode. That being said, it is one of the best possible guardrails against bad data, so consider leaving it to the default value.


## Removed from recent(ish) versions

decontam_use_CDC_varpipe_ref: "If true, use CDC varpipe decontamination reference. If false, use CRyPTIC decontamination reference."  
--> CDC uses their own version of clockwork's decontamination reference, which I call "CDC varpipe" since I pulled it from the varpipe repo. This is currently a null op in myco_raw as it'd require I maintain double the number of Docker images, and it doesn't delineate between human vs NTM vs other forms of contamination (which the task currently requires for some outputs). If there is a demand for CDC varpipe I can make this an option again, but I nevertheless gently recommend against using it due to unclear provenance and contents.

subsample_seed: "Seed used for subsampling with seqtk"  
--> To reduce the number of inputs, this is now hardcoded to 1965.  

tbprofiler_on_bam: "If true, run TBProfiler on BAMs"  
---> No longer supported, we now use Theiagen's fork of TBProfiler for both myco_raw and myco_sra, unless `just_like_2024`


  
# Task-level inputs 
Many of myco's tasks have exposed runtime attributes. These don't show up on the workflow page, but do show up in Terra, and can be accessed by any WDL executor by referencing them in the input JSON. Most of them are for configuring runtime attributes: https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/


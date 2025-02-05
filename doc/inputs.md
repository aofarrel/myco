> [!IMPORTANT]  
> You are currently on an outdated branch of myco that exists solely for reproducing published results. It is HIGHLY recommended you use [a more recent version](https://github.com/aofarrel/myco) in order to take advantage of new updates to clockwork, TBProfiler, and other dependencies.

### A quick introduction to input variables in WDL
There are two kinds of input user-accessible variables in WDL: Workflow-level inputs and task-level inputs. If you are using Terra, you probably don't need to know anything about the difference between them except that task-level inputs get alphabetically sorted in Terra's UI below workflow-level inputs. 

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
It is recommend that you also keep an eye on the total size of your FASTQs. Individual files over subsample_cutoff (default 450 MB, -1 disables this check) will be downsampled, but keep an eye on the cumulative size of samples. For example, a sample like SAMEA968096 has 12 run accessions associated with it. Individually, none of these run accessions' FASTQs are over 1 GB in size, but the sum total of these FASTQs could quickly fill up your disk space. (You probably should not be using SAMEA968096 anyway because it is in sample group, which can cause other issues.)

myco_simple expects that the FASTQs you are putting into have already been cleaned and decontaminated, but this isn't a hard requirement. What is a hard requirement is that you have precisely one forward read and one reverse read per sample -- if you have multi-lane samples across various fastqs, they will need to be merged first.


## Timeouts  
Earlier versions of myco_raw and myco_sra had manual timeouts for the decontamination and variant calling tasks, as timing out was a good way to prevent runaway cloud costs. This is now handled by [guardrail mode](./doc/guardrail_mode.md), with the exception of task-level optional variable `timeout_minutes` for myco_sra's pull task.


## Other variables of note
`decontam_use_CDC_varpipe_ref`: Don't set this to True unless you know what you're doing. This switches the decontamination reference generated from v0.11.3 of clockwork into a decontamination reference [built from the CDC NCHHSTP-DTBE-Varpipe-WGS repo](https://github.com/CDCgov/NCHHSTP-DTBE-Varpipe-WGS/blob/1227ab394a26be0c52a0c9b90c349198681f4f4e/tools/clockwork-0.11.3/OUT/build_references.sh). The history and contents of this reference are not well known, but it is definitely smaller than the typical clockwork v0.11.3 reference. Additionally, its decontamination TSV doesn't specify contaminant types, so instead of knowing (for example) 500 reads were removed for aligning to NTM and 50 for human, you'd just be told 550 reads were contaminants.

*If you want to use v0.12.3+ of clockwork's decontamination reference, which uses CHM13 instead of hg38 for human, please use a newer version of myco. Remember, you're on the archived version here!*

`guardrail_mode`: See (these docs)[./guardrail_mode.md] for more info.

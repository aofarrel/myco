# Quality Control / Sample filtering
Because there is so much QC going on in this pipeline, it can be helpful to divide each QC step into two catagories:
* Site-specific QC: Forms of QC that will remove particular **parts** of a sample. 
* Sample QC: Forms of QC that will remove **an entire sample** from the downstream analysis.

Due to skipping the decontamination step, many filters (site-specific and whole-sample) do not apply to myco_simple.

## Site-specific filtering
A lot of site-specific filtering happens in the decontamination task:

* If over `subsample_cutoff` megabytes, sample is downsampled with seqtk to `subsample_reads`
    * deterministic due to hardcoded seed
    * in myco_sra this happens in the pull task instead of the decontamination task to reduce on file localization costs
* If a reads average quality score is < `fastp_avg_qual`, fastp will remove that read/pair
* Being considered contamination by clockwork

Sites can also be filtered in the variant caller, as clockwork uses minos to adjudicate between a samtools-generated VCF and a cortex-generated VCF to generate a final VCF, as well as in the VCF-to-diff task when they are below `site_min_depth` or located in `call_as_reference_bedfile`. **For people who are more interested in VCFs than MAPLE diff files: Be aware this means your VCF is skipping some filtering steps.**

## Sample filtering

### Typical get-rid-of-the-whole-sample filters (myco_sra and myco_raw)
TBProfiler will filter out samples if the samples average (mean) depth is below `sample_min_avg_depth` or below `sample_min_pct_mapped`

With the exception of `sample_min_q30` samples do not get filtered out in decontamination nor the variant caller unless the task errors out with a handled exception, or a `guardrail_mode` timeout triggers.

A sample can be filtered in VCF-to-diff if more than `sample_max_pct_masked` percent of it is explicitly masked (usually due to low coverage). This means it is informed by site-specific filtering. Note that sites within `call_as_reference_bedfile` do not count as masked with regard to `sample_max_pct_masked` because those sites are considered reference, not an explicit mask.


### FASTQ download (myco_sra only, upstream of all other tasks)
| status code (not reported to user) | situation                                                    | togglable?         | can crash pipeline?         |
|---------------------------|-----------------------------------------------------------------------|--------------------|-----------------------------|
| SRA_BAD_BIOSAMPLE_ID      | BioSample accession appears to be invalid                             | no                 | no                          |
| SRA_FAIL_TO_DOWNLOAD_ALL  | ALL of a BioSample's run accessions fail prefetch and/or fasterq-dump | no                 | `pull.fail_on_invalid`=true |
| SRA_FAIL_TO_DOWNLOAD_SOME | ≥1 of a BioSample's run accessions fail prefetch and/or fasterq-dump  | yes (default: off) | `pull.fail_on_invalid`=true |
| SRA_ONE_FASTQ_ALL         | ALL of a BioSample's run accessions have only one fastq               | no                 | `pull.fail_on_invalid`=true |
| SRA_ONE_FASTQ_SOME        | ≥1 of a BioSample's run accessions have only one fastq                | yes (default: off) | `pull.fail_on_invalid`=true |

Notes: 
* myco_sra does not support sample-level status code outputs; they are defined only in documentation for ease of writing
* SRA_FAIL_TO_DOWNLOAD_ALL usually means that data is corrupt, but it could also mean your network is having issues or SRA is having an outage. If your data is also on ENA, you can try ENABrowserTools or [my WDlization of it](https://github.com/aofarrel/enaBrowserTools-wdl).




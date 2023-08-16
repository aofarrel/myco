# Quality Control / Sample filtering
Not all data is "good" data. Maybe the data is too heavily contaminated, maybe the sequencing machine had trouble with it, maybe its format is incompatiable with our downstream tools, or maybe it's actually a boar from the French countryside that somehow ended up tagged as Mycobacterium tuberculosis. These sorts of things (yes, even the boar one) can and have happened while developing this pipeline, leading to the development of several metrics to help us quickly find and filter out bad data without crashing the entire pipeline. This pipeline was developed further with input from CDPH, who wish to implement QC standards similar to those of the CDC.

Beause there is so much QC going on in this pipeline, it can be helpful to divide each QC step into two catagories:
* Site-specific QC: Forms of QC that will remove particular **parts** of a sample. 
* Sample QC: Forms of QC that will remove **an entire sample** from the downstream analysis.

## Site-specific filtering

### subsampling
If a FASTQ is above `subsample_cutoff` MB, it will get downsampled by seqtk. `subsample_cutoff` is turned off (set to -1) by default in myco_raw, and set to 450 by default in myco_sra. myco_cleaned does not support subsampling.

### decontamination
An entire WDL task of myco (except myco_cleaned) is dedicated just to decontaminating reads. The decontamination workflow starts with `clockwork map_reads` to map to a decontamination reference, and then uses `clockwork rm_contam` to generate decontaminated FASTQs. It is worth noting that how long a sample spends in this decontamination step roughly correlates with how much contamination is in it, but input file size is also a factor. If you're seeing a batch of samples that are roughly the same size (or subject to default downsampling settings) as typical, but take unusually long to decontaminate, that batch of samples might be considered suspect.

### earlyQC (aka TBfastProfiler)
EarlyQC merges TBProfiler (in fastq-input-mode) and fastp into one WDL step which will run unless `early_qc_skip_entirely` is true or you are running myco_cleaned. TBProfiler does no site-specific filtering of its own, but if `early_qc_skip_trimming` is false, fastp will further clean your FASTQs as a form of site-specific filtering. 

#### removing low-quality read pairs
`early_qc_trim_qual_below` is piped into fastp's `average_qual`. If one read's average quality score is < `average_qual`, then that read/pair is discarded. You can disable this by setting `early_qc_trim_qual_below` to 0 or `early_qc_skip_trimming` to true.

### variant calling
The variant caller used by all forms of myco uses clockwork, which itself leverages minos. minos will generate VCFs using two different methods, then compare the two of them, then output a final ajudicated VCF.

### VCF to diff
Lily Karim's VCF to diff script will mask any sites that are below `diff_min_site_coverage` (default: 10) from appearing in the final diff file output. It does not affect the VCF.

## Sample filtering
Generally speaking, these filters apply to myco_sra and myco_raw. The only sample-level filtering myco_cleaned supports are TREE_TOO_MANY_LOW_COVERAGE_SITES, VARIANT_CALLING_KILLED, and VARIANT_CALLING_TIMEOUT.

### FASTQ download (myco_sra only)
| status code               | situation                                                             | togglable?         | can crash pipeline?         |
|---------------------------|-----------------------------------------------------------------------|--------------------|-----------------------------|
| SRA_BAD_BIOSAMPLE_ID      | BioSample accession appears to be invalid                             | no                 | no                          |
| SRA_FAIL_TO_DOWNLOAD_ALL  | ALL of a BioSample's run accessions fail prefetch and/or fasterq-dump | no                 | `pull.fail_on_invalid`=true |
| SRA_FAIL_TO_DOWNLOAD_SOME | ≥1 of a BioSample's run accessions fail prefetch and/or fasterq-dump  | yes (default: off) | `pull.fail_on_invalid`=true |
| SRA_ONE_FASTQ_ALL         | ALL of a BioSample's run accessions have only one fastq               | no                 | `pull.fail_on_invalid`=true |
| SRA_ONE_FASTQ_SOME        | ≥1 of a BioSample's run accessions have only one fastq                | yes (default: off) | `pull.fail_on_invalid`=true |

Notes: 
* myco_sra does not support sample-level status code outputs; they are defined only in documentation for ease of writing
* SRA_FAIL_TO_DOWNLOAD_ALL usually means that data is corrupt, but it could also mean your network is having issues or SRA is having an outage. If your data is also on ENA, you can try ENABrowserTools or [my WDlization of it](https://github.com/aofarrel/enaBrowserTools-wdl).
* These timers apply to the same WDL task but are for different processes within that task -- `timeout_decontam_part1` is 20 and `timeout_decontam_part1` is 15, and a sample spends 19 minutes mapping plus another 14 minutes finishing the decontamination process, it will *not* be filtered out.

### decontamination
Entire samples do not get filtered out here unless the decontamination task errors out, or you have timeouts -- specifically `timeout_decontam_part1` and `timeout_decontam_part2` -- set to a nonzero value. The reason for timeouts filtering out samples is that a sample taking a long time is itself a sign that the sample is heavily contaminated, and a heavily decontaminated sample is more likely to have too many sites removed for variant calling to work properly, which is useful if processing tens of thousands of samples from SRA of varying degrees of quality. It is, however, a lot fuzzier than most other forms of QC in this pipeline, so timeouts are turned off (set to 0) by default for myco_raw. For more information on the circumstances that can cause the decontamination task to error out, please see [status_codes.md](./status_codes.md).

### earlyQC 
If more than `early_qc_minimum_q30` (as float where 0.5=50%) of your decontaminated FASTQs's calls are below Q30, and if `early_qc_skip_qc` is false, and if `early_qc_skip_qc` is also false, the sample will be removed with status `EARLYQC_TOO_MANY_BELOW_Q30`. This is independent of fastp's site-specific filtering (eg, `early_qc_skip_trimming`, and `early_qc_trim_qual_below`).

### variant calling
As with decontamination, entire samples do not get filtered out here unless the variant caller has an error or times out.

### covstats
TODO

### VCF to diff
The site-specific filtering of VCF-to-diff informs the sample-level filtering. If too many of a sample's sites are masked for having coverage below `diff_min_site_coverage`, the entire sample will be throw out instead. "Too many" is defined as float `diff_max_pct_low_coverage`, where 0.5=50%. You can effecitvely turn off this filter by setting `diff_max_pct_low_coverage` to 1.01 or higher.

### Tree Nine
myco_cleaned does not support VCF-to-diff's ability to filter out samples based on having too many low coverages sites, so it technically does that filtering in Tree Nine instead. It's the same filter, just done in a different task. To avoid redundancy, myco_raw and myco_sra do not allow the user to attempt to do this filtering in Tree Nine, since for them it's already being done earlier.




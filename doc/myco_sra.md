# myco_sra
myco_sra is the [SRA](https://www.ncbi.nlm.nih.gov/sra) version of myco. Use this version of myco if you want to analysze fastqs from SRA. This is powered by [SRANWRP](https://github.com/aofarrel/SRANWRP), most notably the [pull-FASTQs-from-biosample](https://dockstore.org/workflows/github.com/aofarrel/SRANWRP/pull_FASTQs_from_SRA_by_biosample:main?tab=info) workflow. (For the sake of simplicity this readme calls the part of myco_sra that downloads from SRA "SRANWRP", even thought SRANWRP contains a few additional utility functions and workflows.)

Input: biosample_accessions -- a single text File which contains BioSample accessions, one BioSample accession per line, ex:
```
SAMEA10030079
SAMEA10030285
SAMEA10030321
SAMEA10030646
SAMEA104390589
SAMEA110024138
SAMEA111556114
```

SAME, SAMN, SRS, ERS, and numeric BioSample accessions are all supported, as well as any combination of these formats. **Run accessions (SRR, ERR, DRR) are not supported, [but you can use this workflow to convert your run accessions to BioSample accessions](https://dockstore.org/workflows/github.com/aofarrel/SRANWRP/get_biosample_accessions_from_run_accessions:main?tab=info).**

All other inputs are documented here: [doc/inputs.md](./doc/inputs.md)

## How does the fastq downloading part differ from similar workflows?
There are several existing workflows which can pull from SRA, such as [SRA Fetch](https://dockstore.org/workflows/github.com/theiagen/terra_utilities/SRA_Fetch:v1.4.1?tab=info) and [DownloadFromSRA](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/DownloadFromSRA:kvg_update_downloaders?tab=info). If you need your reads downloaded with no processing or saved to a specific GCS directory, these workflows might be better suited to your needs than SRANWRP, because SRANWRP assumes you only want paired-end Illumina data while also making no assumptions about you actually giving it that.

SRANWRP also supports downloading fastqs per BioSamples accession rather than per run accession. In myco_sra this helps maintain a sense of which reads belong to which sample (clockwork operates per-sample), but it can also be helpful when dealing with samples with very large number of run accessions.

## How does myco_sra handle "weird" data?/What checks does it perform?
The "ideal" scenario is that each BioSample has some number of run accessions, and each run accession returns one pair of Illumina-created fastq files. But sometimes life isn't so simple, so SRANWRP handles:

* BioSamples which contain run accessions created with different sequencing technologies (ex: SAMN03257097 has two Illumina run accessions and two PacBio run accessions)
* Run accessions that return three fastq files
* Run accessions that return one fastq file
* Run accessions that return more than one pair of fastqs
* Run accessions that return large fastqs -- by default, fastqs larger than 450 megabytes will be randomly subsampled down to 1,000,000 random reads
* BioSamples with any combination of the above exist -- ex: A single BioSample with one run accession returning one pair of fastqs, one run accession returning three fastqs, one run accession actually being a PacBio run, and one run accession returning two pairs of 2 GB fastqs

## Are there any SRA accessions known to break SRANWRP?
Accessions belonging to "sample groups" are not supported, as it isn't very clear which run accessions correlate to which sample accessions. Known "sample group" accessions can be found [in SRANWPR's denylists](https://github.com/aofarrel/SRANWRP/tree/main/inputs/denylists). Non-TB accessions are also not supported.

Although none that aren't also in sample groups have been found yet, in theory, samples with a very large number of run accessions could be problematic on a GCP backend. GCP backends require you request a certain amount of disk size before runtime, which can get dicey when what you want to do at runtime is download an unknown number of files of unknown size. As such, SRANWRP requests more disk size than you *probably* need, but there could come a time that guess doesn't prove to be enough.

The only other accessions known to fail SRANWRP are ones which break prefetch or fasterq-dump. Usually this means the fastqs themselves are invalid (ex: SAMEA1877221, SAMEA2609926, SAMEA2609935) or were not fully processed by NCBI (ex: ERR760606, although SAMEA3231653 has two other run accessions which do work fine).

## What if an SRA accession passes SRANWRP, but isn't good enough for the rest of the pipeline?
We can't always tell that data isn't up to our standards until later down the pipeline. myco (including myco_sra) will filter out samples which:
* are too heavily contaminated to complete decontamination in a timely manner<sup>†</sup>
* take too long in the variant caller<sup>†</sup>
* have low overall coverage
<sup>†</sup>This sort of filtering can be disabled by setting the timeout optional variables to 0 -- but be aware that GCP will kill any VM that is still alive after about a week, so if you're on Terra, your samples need to process faster than that!

There is a small number of BioSamples known to return fastqs which, after decontamination, will cause a runtime error in the variant caller. As of version 3.0.1 of myco, this runtime error is handled by throwing out the sample rather than stopping the entire pipeline (unless you set crash_on_error to true). Throwing out samples like this was decided as the default behavior as a quick analysis indicates the runtime error only happens when clockwork's variant calling removes >95% of the fastq during a Trimmomatic run (ie, the sample probably shouldn't be used anyway). The list of these samples can be found [in SRANWPR's denylists](https://github.com/aofarrel/SRANWRP/tree/main/inputs/denylists).
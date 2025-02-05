# FAQs

## What's the use case of each version of myco?
Short answer: If you are are an LHJ or anyone else who has non-SRA TB samples lying around, use myco_raw. If your samples are on SRA, use myco_sra. If your samples should undergo as little processing (eg: cleaning, decontamination) as possible and you're okay with having less features, use myco_simple.

Longer answer:
* pairs of FASTQs which have been decontaminated and merged such that each sample has precisely two FASTQs associated with it: **myco_simple** 
 * pairs of FASTQs which have yet to be decontaminated or merged: 
     * if each sample has its FASTQs in a single array: **myco_raw** 
     * if each sample has its forward FASTQs in one array and reverse FASTQs in another array: [Decontam_And_Combine_One_Samples_Fastqs](https://dockstore.org/workflows/github.com/aofarrel/clockwork-wdl/Decontam_And_Combine_One_Samples_Fastqs), then **myco_simple** 
     * a list of SRA BioSamples whose FASTQs you'd like to use **myco_sra** 
     * a list of SRA run accessions (ERR, SRR, DRR) whose FASTQs you'd like to use: [convert them to BioSamples](https://dockstore.org/workflows/github.com/aofarrel/SRANWRP/get_biosample_accessions_from_run_accessions:main?tab=info), then **myco_sra**)   


## When a sample is filtered out, does the entire pipeline/workflow/Terra run crash?
On default settings, myco_raw and myco_sra will attempt to throw out bad samples silently, without causing the workflow to exit with an error. This is great for people who want to run as many samples as possible at once, because one bad sample will not cause everything else to halt. It's also useful for those who use Terra data tables, as intermediate outputs and a status code can be saved to the Terra data table, making it possible to keep some intermediate outputs and see at a glance which samples are problematic.

## How do I know where a sample got filtered out?
If you are using myco_raw such that each instance of myco_raw receives only one sample -- if you use Terra data tables where each row represents one sample, this means you! -- you can get info on where your sample got dropped from the status code output. See status_codes.md for more info and how to interpret status codes.

If you are running any of the pipelines such that each instance of the workflow recieves more than one sample, you will need to interpret your output files to see where samples got dropped. For instance, a sample that has a VCF but no diff must have failed in the vcf-to-diff task.

See qc_and_filtering.md for more information on filtering out samples.

## Can myco filter out bad data within a sample without throwing out the whole sample?
Yes. See qc_and_filtering.md for more information.

## Why are there so many places where samples get dropped?
Myco was originally designed with the goal of analyzing as much TB SRA data as we could relatively quickly and cheaply. About 94% of MTBC BioSamples which are tagged as containing paired Illumina data pass myco_sra default's filters. An additional 3% of MTBC BioSamples are either completely invalid (fails fasterq-dump due to a database migration error that happended years ago, length of quality score doesn't match length of nucleotide strings, etc) or can't be used by clockwork (not actually paired Illumina reads, etc). The remainder of what gets filtered out seems to be extremely heavily contaminated and/or low overall coverage. As such, we're reasonably confident that myco_sra's default filters are a good tradeoff for the realities of SRA's data. Of course, you might want to cast a wider net than we did, or you might not be using SRA data at all, so all of these filters except for this-will-100%-break-the-pipeline-if-you-let-this-through ones can be configured.

## What if I want to use a different reference genome?
This isn't officially supported due to TBProfiler, UShER, and clockwork each needing specific reference genomes:
* TBProfiler's reference genome must be *exactly* the same as the one you called variants upon, if you're running TBProfiler on bams
* UShER has limitations on how long a chromosome's name can be
* clockwork's decontamination is designed with their specific decontamination reference in mind
* covstats also has a reference genome in mind
That being said, old versions of myco used clockwork reference prepare to prepare the TB genome, and you could hack that passing-reference-genomes-around functionality to use your own custom reference genomes if you're confident. The latest version of myco_raw and myco_sra that used clockwork reference prepare was [4.1.3](https://github.com/aofarrel/myco/releases/tag/4.1.3), so that's a good place to start. The latest version of myco that required reference genomes for TBProfiler and/or UShER was [4.2.0](https://github.com/aofarrel/myco/releases/tag/4.2.0). This version of myco does not have sample-level status codes.

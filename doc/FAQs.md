# FAQs

<!--### How do I install and run TB-D?
If you intend on using TB-D on the Terra website, please see [Get Started on Terra](./get_started_Terra.md).
For all other users, please see [this webpage](.get_started_nonTerra.md) instead.-->

## Basic questions

### What sequencing technologies are supported?
This pipeline relies on clockwork, which only supports PE Illumina. 

### What's the use case of each version of myco?
Short answer: If you are are an LHJ or anyone else who has non-SRA TB samples lying around, use myco_raw. If your samples are on SRA, use myco_sra. If your samples should undergo as little processing (eg: cleaning, decontamination) as possible and you're okay with having less features, use myco_simple.

Longer answer:
* pairs of FASTQs which have been decontaminated and merged such that each sample has precisely two FASTQs associated with it: **myco_simple** 
  * if these are in Terra data table format, you may want to use the **wrapper_example** 
 * pairs of FASTQs which have yet to be decontaminated or merged: 
     * if each sample has its FASTQs in a single array: **myco_raw** 
     * if each sample has its forward FASTQs in one array and reverse FASTQs in another array: [Decontam_And_Combine_One_Samples_Fastqs](https://dockstore.org/workflows/github.com/aofarrel/clockwork-wdl/Decontam_And_Combine_One_Samples_Fastqs), then **myco_simple** 
     * a list of SRA BioSamples whose FASTQs you'd like to use **myco_sra** 
     * a list of SRA run accessions (ERR, SRR, DRR) whose FASTQs you'd like to use: [convert them to BioSamples](https://dockstore.org/workflows/github.com/aofarrel/SRANWRP/get_biosample_accessions_from_run_accessions:main?tab=info), then **myco_sra**)   


## Inputs
*More info: [inputs.md](./inputs.md)*  

### How should I enter my fastqs?
myco_sra takes in a text file of NCBI BioSample accessions, one per line. BioSample accessions may be in SAMN/SAME/SAMD format, SRS/ERS/DRS format, numerical ID, or any combination thereof. In accordance with the NCBI data model and clockwork's recommendations, a single BioSample is treated as a single sample, regardless of the number of run accessions (SRR/ERR/DRR) are within it, but be aware some submitters do not follow this data model. For example, some people upload 100 run accessions, representing 100 samples, under a single BioSample.

myco_raw and myco_simple take in files directly as a nested array, where every inner array is one sample's fastq files. It is expected that every valid sample has 2, 4, or some other even number of fastq files.

### Does myco support multi-lane/multi-run samples?
Yes. myco_sra and myco_raw will automatically concatenate multi-lane samples into a R1 and an R2 file. myco_simple on the other hand needs to be given precisely one R1 FASTQ file and precisely one R2 FASTQ file.

For example, SAMN02584599 contains both SRR1166330 and SRR1169013. When pulling via myco_sra, you will end up with SRR1166330_1.fastq, SRR1166330_2.fastq, SRR1169013_1.fastq, and SRR1169013_2.fastq in the working directory. These get renamed to SAMN02584599_SRR1166330_1.fastq, SAMN02584599_SRR1166330_2.fastq, SAMN02584599_SRR1169013_1.fastq, and SAMN02584599_SRR1169013_2.fastq in order to keep the name of the sample in the filenames. These files are then passed to the decontamination task, and assuming they pass, your cleaned FASTQs will be called SAMN02584599_1.fastq and SAMN02584599_2.fastq and then carry on to the variant caller. The output of the variant caller will be SAMN02584599.vcf, from which you get SAMN02584599.diff, which will appear on your tree as SAMN02584599.

### Will myco_sra crash if I accidentally pass in accessions of non-Illumina / corrupted samples?
Unless the samples are bad in a way I failed to account for, no. myco_sra includes a wrapper of an intentionally outdated version of fasterq-dump that only supports Illumina files. The wrapper also does its best to handle barcode files, multi-lane samples, mixed-sequencing-technology-files, and all other edge cases I could handle. Samples that cannot be pulled will be silently skipped. Do note that myco_sra cannot detect if a sample is in a sample pool, has no idea if a sample is *in silico* or a resequence, and cannot predict file size prior to pulling; see notes in [inputs.md](./inputs.md) for more info.

### How should I name my files to make sure myco_raw/myco_simple knows which file is R1 and which one is R2?
If all of your samples have precisely one R1 file and precisely one R2 file, we recommend SAMPLE_R1.fastq and SAMPLE_R2.fastq format. If any of your inputs are multi-lane samples that have not been concatenated, we recommend SAMPLE_RUN_1.fastq and SAMPLE_RUN_2.fastq. These are not strict requirements -- other iterations (`_1/_2`, `.fq`, etc) may work, but do a small-scale test first (or just rename your samples).

## Filters

### Can myco filter out bad data within a sample without throwing out the whole sample?
Yes, this is done via fastp cleaning, decontamination, and masking of low-coverage sites.

### When a sample is filtered out, does the entire pipeline/workflow/Terra run crash?
On default settings, myco will attempt to throw out bad samples silently, without causing the workflow to exit with an error. This is great for people who want to run as many samples as possible at once, because one bad sample will not cause everything else to halt. It's also useful for those who use Terra data tables, as intermediate outputs and a status code can be saved to the Terra data table, making it possible to keep some intermediate outputs and see at a glance which samples are problematic.

To cause most handled exceptions (but not mere QC failures) to crash, set `decontam_each_sample.crash_loudly`, `variant_calling.crash_on_error`, and `variant_calling.crash_on_timeout` to true (all default to false). Handled exceptions include scenarios such as the variant caller outputting an empty VCF or running out of memory.

### What kind of samples will myco filter out?
myco will filter out samples that fail its QC standards. Additionally, when downloading from NCBI SRA, TB-D_sra will filter out:
* [Index/Barcode sequences](https://www.biostars.org/p/390726/#390738) that are sequenced seperately from R1/R2
  * If a sample downloads as three FASTQ files, then that index/barcode file will be thrown out and the remaining R1/R2 pair will be kept
  * If a sample downloads as one, five, seven, nine, etc FASTQ files, the whole sample will be thrown out
* Non-Illumina samples (PacBio, etc)
* SE Illumina samples
* Corrupt files, such as when the FASTQ quality score isn't the same length as the nucleotide string

Please be aware that myco cannot automatically detect:
* *in silico* synthetic data such as SAMN18146425
* samples within "sample pools" such as SAMEA968074
* samples that contain data from multiple biological samples that were incorrectly submitted as single-sample

...but you might be able to spot such oddities on a phylogenetic tree.

### When a sample is filtered out, does the entire pipeline/workflow/Terra run crash?
On default settings, myco_raw and myco_sra will attempt to throw out bad samples silently, without causing the workflow to exit with an error. This is great for people who want to run as many samples as possible at once, because one bad sample will not cause everything else to halt. It's also useful for those who use Terra data tables, as intermediate outputs and a status code can be saved to the Terra data table, making it possible to keep some intermediate outputs and see at a glance which samples are problematic.

You can set `variantcalling_crash_on_error` to true (default: false) to make the pipeline crash if any error is encountered in the variant caller. "Any error" includes the variant caller running out of memory, timing out, getting an unknown error, or (most relevant to data QC) failing to actually call enough variants due to the sample being too small/corrupt/contaminated. In other words, the variant caller task getting an error is a strong sign that the input data is of questionable quality, so you may want to crash there if you are unsure about the quality of your samples and would prefer to cut your losses early if one of them proves suspect.

### How do I know where a sample got filtered out?
If you are using myco_raw such that each instance of myco_raw receives only one sample -- if you use Terra data tables where each row represents one sample, this means you! -- you can get info on where your sample got dropped from the status code output. See status_codes.md for more info and how to interpret status codes.

If you are running any of the pipelines such that each instance of the workflow recieves more than one sample, you will need to interpret your output files to see where samples got dropped. For instance, a sample that has a VCF but no diff must have failed in the vcf-to-diff task.

See qc_and_filtering.md for more information on filtering out samples.

### Why are there so many places where samples get dropped?
myco was originally designed with the goal of analyzing as much TB SRA data as we could relatively quickly and cheaply. About 94% of MTBC BioSamples which are tagged as containing paired Illumina data pass myco_sra default (for that version) filters. An additional 3% of MTBC BioSamples are either completely invalid (fails fasterq-dump due to a database migration error that happended years ago, length of quality score doesn't match length of nucleotide strings, etc) or can't be used by clockwork (not actually paired Illumina reads, etc). The remainder of what gets filtered out seems to be extremely heavily contaminated and/or low overall coverage. As such, we're reasonably confident that myco_sra's default filters are a good tradeoff for the realities of SRA's data. Of course, you might want to cast a wider net than we did, or you might not be using SRA data at all, so all of these filters except for this-will-100%-break-the-pipeline-if-you-let-this-through ones can be configured.

## Atypical use cases

### I want to replicate your published results as closely as possible. Which version of TB-D should I use?
Several new features, user-friendly options, and critical bugfixes have been added to the pipeline since it was used to generate data, so older versions of the pipeline aren't supported, including the "publication reproducibility" branch. However, the general logic flow of the pipeline has remained almost identical, so if you set `just_like_2024` to True to change the few "important" differences back to how they were before:
* Older version (0.11.3) of the decontamination reference
* Older version (0.11.3) of Clockwork
* TBProfiler is run on BAM files instead of FQs
* Some QC and subsampling adjustments

### What if I want to use a different reference genome?
This isn't officially supported due to TBProfiler, UShER, and clockwork each needing specific reference genomes:
* TBProfiler's reference genome must be *exactly* the same as the one you called variants upon, if you're running TBProfiler on bams
* UShER (in Tree Nine, which is often run after myco) has limitations on how long a chromosome's name can be
* clockwork's decontamination is designed with their specific decontamination reference in mind
* covstats also has a reference genome in mind
That being said, old versions of myco used clockwork reference prepare to prepare the TB genome, and you could hack that passing-reference-genomes-around functionality to use your own custom reference genomes if you're confident. The latest version of myco_raw and myco_sra that used clockwork reference prepare was [4.1.3](https://github.com/aofarrel/myco/releases/tag/4.1.3), so that's a good place to start. The latest version of myco that required reference genomes for TBProfiler and/or UShER was [4.2.0](https://github.com/aofarrel/myco/releases/tag/4.2.0). This version of myco does not have sample-level status codes.

## Troubleshooting

### My tasks are getting cancelled due to lack of resources (sigkill, return code 137, out-of-memory, etc) and I am not running on Terra
Your system resources may not be may not be enough to run the pipeline as intended. Usually the "limiting reagent" is RAM -- we recommend a minimum of 32 GB, although 16 GB will generally work, especially if you enable downsampling. Be aware that Docker Desktop by default will throttle how much memory a Docker container will have, so you might be giving tasks less than you expect!

It's worth noting that on non-cloud backends, Cromwell is not good at allocating resources for concurrent tasks, so a workflow that runs fine on miniwdl might crash on Cromwell. You can sort of sidestep this by adjusting Cromwell's configuration file so concurrent-job-limit is 1, but it's probably better to just switch to miniwdl, which (usually) properly handles concurrent tasks automatically and without crashing.

More info:
* [Cromwell/Dockstore CLI only] [Troubleshooting the most common resource issues](https://docs.dockstore.org/en/stable/advanced-topics/dockstore-cli/dockstore-cli-faq.html#cromwell-docker-lockup) -- this document is in the context of the Dockstore CLI, which runs Cromwell under the hood, but it also applies to Cromwell itself
* [Resource constraints on Docker Engine](https://docs.docker.com/engine/containers/resource_constraints/)
* [Resource constraints on Docker Desktop](https://docs.docker.com/desktop/settings-and-maintenance/settings/#advanced)

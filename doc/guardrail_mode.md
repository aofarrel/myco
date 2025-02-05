> [!IMPORTANT]  
> You are currently on an outdated branch of myco that exists solely for reproducing published results. It is HIGHLY recommended you use [a more recent version](https://github.com/aofarrel/myco) in order to take advantage of new updates to clockwork, TBProfiler, and other dependencies.

# Guardrail Mode
Guardrail mode is designed to filter out problematic samples that may cause the pipeline to take several hours longer than expected and/or crash. This is done by setting several "safeguards" which will discard a sample if it's really, really bad. To make it easier for users, safeguards are not directly adjustable by the user -- instead, the user simply turns guardrail mode on or off.

Guardrail mode is **not** a replacement for QC. It is arguably a type of QC filter, yes, but it is an extremely lenient one that is designed to remove only samples which cannot be meaningfully analyzed.


## Rationale
When working with data of unknown quality, it can be helpful to quickly remove samples that are likely low-quality. While developing myco on SRA data, we noticed that if a given sample took an unusually long time in the decontamination or variant calling step, they were likely to end up filtered out by the final quality control steps of the pipeline. This is especially true of the decontamination step -- the more contamination a sample has, the more that step has to do. This heuristic was defined on the default runtime attributes and using Terra as a backend, so straying from those defaults -- including changing from SDDs to HDDs -- is likely to make the default timeout values less useful.


## Filters
When guardrail mode is active, if any of these are true about a sample, the sample will be removed.
* TBProfiler thinks the median coverage is less than 3x or that more than 10% of the data can't map to H37Rv
* Mapping reads to the decontamination reference takes more than 300 minutes
* Decontaminating takes more than 600 minutes
* Variant calling takes more 600 minutes
* fastp determines that a less than 20% of the sample has a q score above 30, as measured before fastp cleaning and before decontam

Additionally, specific to myco_sra, if any input fastq is larger than `subsample_cutoff` GB, it will be heavily downsampled. You can turn off this behavior by setting `subsample_cutoff` to -1.
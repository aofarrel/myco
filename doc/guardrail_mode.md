# Guardrail Mode
Guardrail mode is designed to filter out problematic samples that may cause the pipeline to take several hours longer than expected and/or crash. This is done by setting several "safeguards" which will discard a sample if it's really, really bad. To make it easier for users, safeguards are not directly adjustable by the user -- instead, the user simply turns guardrail mode on or off.

Guardrail mode is **not** a replacement for QC. It is arguably a type of QC filter, yes, but it is an extremely lenient one that is designed to remove only samples which cannot be meaningfully analyzed.


## Filters
When guardrail mode is active, if any of these are true about a sample, the sample will be removed.
* TBProfiler thinks the median coverage is less than 3x or that more than 10% of the data can't map to H37Rv
* Mapping reads to the decontamination reference takes more than 300 minutes
* Decontaminating takes more than 600 minutes
* Variant calling takes more 600 minutes
* fastp determines that a less than 20% of the sample has a q score above 30, as measured before fastp cleaning and before decontam

Additionally, if any input fastq is larger than `subsample_cutoff` GB, it will be heavily downsampled.


## Are there really samples out there that would fail these filters?
Yes. All filters implemented by Guardrail mode represent actual issues I had with SRA data. The Q30 filter is the only exception -- the verison of the pipeline that ran on SRA data did not have fastp, so that one is more of a hypothetical.
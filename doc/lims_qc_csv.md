## QC CSV Fields
sample
    Sample ID
genome_pct_coverage
    Percentage of the genome covered by the sample (earlyQC)
mean_coverage
    Mean coverage (covstats)
median_coverage
    Median coverage (earlyQC, specifically TBProfiler)
pct_above_q30
    Percent of reads above q30 (earlyQC, specifically fastp)
reads_is_contam
    Reads considered contaminated (and removed) by the decontaminator
reads_reference
    Reads considered mapping to reference by the decontaminator
reads_unmapped
    Reads considered unmapped by the decontaminator (NOT the variant caller!)
status
    The status of this sample
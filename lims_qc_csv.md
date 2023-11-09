
## implemented
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
    Reads considered unmapped by the decontaminator
status
    The status of this sample

## ideas
warnings
    Non-fatal warnings for this sample
out_pct_pass_fastp
    Percent of reads that pass fastp's filters
coverage_cutoff
    Coverage below this is considered "low coverage" wrt in_pct_low_cov
in_pct_low_cov
    Percent of sites (CHECK!!) determined with coverage below coverage_cutoff; these sites get masked when creating diff files

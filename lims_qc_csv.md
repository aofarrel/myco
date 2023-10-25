draft draft draft

sample
    Sample ID
status
    The status of this sample
warnings
    Non-fatal warnings for this sample
    
reads_is_contam
    Reads considered contaminated (and removed) by the decontaminator
reads_reference
    Reads considered mapping to reference by the decontaminator
reads_unmapped
    Reads considered unmapped by the decontaminator

FASTP
in_dupe_rate
out_dupe_rate
in_pct_q30
out_pct_q30
out_pct_pass_fastp
    Percent of reads that pass fastp's filters

median_coverage
    Median coverage as determined by TBProfiler
mean_coverage
    Mean coverage as determined by covstats

coverage_cutoff
    Coverage below this is considered "low coverage" wrt in_pct_low_cov
in_pct_low_cov
    Percent of sites (CHECK!!) determined with coverage below coverage_cutoff; these sites get masked when creating diff files
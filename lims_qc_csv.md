draft draft draft

1. sample: sample ID 
2. eQC_dupe_rate: duplication rate (fastp)
3. eQC_pct_pass: percent of reads that pass fastp's filters
4. eQC_pct_q30
5. dQC_pct_low_cov: percent of variants that are low coverage (vcf-to-diff) (CHECK!)
6. pQC_med_cov: median coverage of the sameple (tbprofiler)
7. cQC_mean_cov: mean coverage of the sample (covstats)
8. input_size: size of the raw fastq inputs
9. decon_size: size of the decontaminated fastqs before earlyQC
10. status: the status/"error code" for this sample
11. warnings: non-fatal warnings for this sample
# status codes
Status codes represent the status of a given sample after it has completed the pipeline, providing information as to whether it completed the pipeline successfully or was filtered out somewhere along the way. They are visible to the user as a workflow-level output only if you are running myco_raw on a Terra data table, where each instance of the myco_raw workflow recieves one sample.

| status code                             | explanation                                                                                                                           | sample can be recovered?       | suggested resolution                                                                                                                                                                                                                                                                                                                                                                                                   |
|-----------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------|--------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| COVSTATS_BAD_MAP_AND_COVERAGE           | Combination of `COVSTATS_LOW_MEAN_COVERAGE` and `COVSTATS_LOW_PCT_MAPPED_TO_REF`                                                      | yes, but the sample is suspect | Your sample will will likely be so heavily masked in VCF-to-diff that it's not worth salvaging, but you can try with `covstats_qc_skip_entirely` = true                                                                                                                                                                                                                                                                |
| COVSTATS_LOW_MEAN_COVERAGE              | Covstats detected that less your sample's mean coverage was less than `covstats_qc_cutoff_coverages`                                  | yes, but the sample is suspect | Your sample will will likely be so heavily masked in VCF-to-diff that it's not worth salvaging, but you can try with `covstats_qc_skip_entirely` = true                                                                                                                                                                                                                                                                |
| COVSTATS_LOW_PCT_MAPPED_TO_REF          | Covstats detected that more than `covstats_qc_cutoff_unmapped` percent of your sample was mapped to the H37Rv reference               | yes, but the sample is suspect | Rerun with a higher value for `covstats_qc_cutoff_unmapped`.                                                                                                                                                                                                                                                                                                                                                           |
| DECONTAMINATION_MAP_READS_TIMEOUT       | The map_reads part of the decontamination process went over `timeout_decontam_part1` minutes                                          | yes, but the sample is suspect | This could be a sign your sample is very heavily contaminated. If you wish to continue attempting to use it, set `timeout_decontam_part1` to 0 and rerun.                                                                                                                                                                                                                                                              |
| DECONTAMINATION_MAP_READS_KILLED        | The map_reads part of the decontamination process was killed (return code 137)                                                        | yes                            | Set the decontamination task's memory runtime attribute to a higher value (default: 16 GB) and rerun.                                                                                                                                                                                                                                                                                                                  |
| DECONTAMINATION_MAP_READS_UNKNOWN_ERROR | The map_reads part of the decontamination process had an unknown error                                                                | no                             | Open an issue on GitHub                                                                                                                                                                                                                                                                                                                                                                                                |
| DECONTAMINATION_RM_CONTAM_KILLED        | The rm_contam part of the decontamination process was killed (return code 137)                                                        | yes                            | Set the decontamination task's memory runtime attribute to a higher value (default: 16 GB) and rerun.                                                                                                                                                                                                                                                                                                                  |
| DECONTAMINATION_RM_CONTAM_UNKNOWN_ERROR | The rm_contam part of the decontamination process had an unknown error                                                                | no                             | Open an issue on GitHub                                                                                                                                                                                                                                                                                                                                                                                                |
| DECONTAMINATION_RM_CONTAM_TIMEOUT       | The rm_contam part of the decontamination process went over `timeout_decontam_part2` minutes                                          | yes, but the sample is suspect | This could be a sign your sample is very heavily contaminated. If you wish to continue attempting to use it, set `timeout_decontam_part2` to 0 and rerun.                                                                                                                                                                                                                                                              |
| EARLYQC_TOO_MANY_BELOW_Q30              | fastp detected `early_qc_cutoff_q30`*100 percent of your FASTQs's calls have a quality score below 30                                 | yes, but the sample is suspect | This could be a sign your sample is very low quality, possibly due issues in sample purification or during sequencing. If you wish to continue attempting to use it, adjust `early_qc_cutoff_q30` to a lower value (default: 0.90)                                                                                                                                                                                     |
| VARIANT_CALLING_KILLED                  | The variant calling task was killed (return code 137)                                                                                 | yes, but the sample is suspect | Set `variantcalling_memory` to a higher value (default: 32 GB) and rerun, but be aware that running out of memory on default settings is quite unusual and may indicate an issue with the data.                                                                                                                                                                                                                        |
| VARIANT_CALLING_TIMEOUT                 | The variant calling task went over `timeout_variant_caller` minutes                                                                   | yes, but the sample is suspect | This could be a sign your sample is very small or very large. If you wish to continue attempting to use it, set `timeout_variant_caller` to 0.                                                                                                                                                                                                                                                                         |
| VARIANT_CALLING_UNKNOWN_ERROR           | The variant calling task returned 1 for unknown reasons                                                                               | no                             | Your FASTQs might be corrupt or almost entirely empty.                                                                                                                                                                                                                                                                                                                                                                 |
| VARIANT_CALLING_UNKNOWN_ERROR_$rc       | The variant calling task returned $rc for unknown reasons                                                                             | no                             | Your FASTQS might be corrupt or almost entirely empty.                                                                                                                                                                                                                                                                                                                                                                 |
| VARIANT_CALLING_ADJUDICATION_FAILURE    | The variant calling task failed, and it appears your sample has enough sites for minimap2 but not Cortex                              | yes, if sample can be bigger   | It appears Cortex cannot find any variants to call. It's possible too much of it was removed during the decontamination step, or there was never much of it in the first place. Check the size of this sample's input FASTQs and compare that to the size of the FASTQs after the decontamination step and earlyQC. You *might* be able to recover this sample by running myco_cleaned on raw, not-downsampled FASTQs. |
| VCF2DIFF_TOO_MANY_LOW_COVERAGE_SITES    | VCF-to-diff task found ≥`diff_max_pct_low_coverage`*100 percent of sample's sites are below `diff_min_site_coverage` coverage         | yes, but the sample is suspect | A diff file can still be generated if `diff_min_coverage_per_site` (default: 10) is set to 0, but note that low coverage sites will not be masked in the resulting diff file.                                                                                                                                                                                                                                                                                                                                                                    

### not visible to user, but defined in documentation
| status code                             | explanation                                                           | sample can be recovered?       | suggested resolution                                   |
|-----------------------------------------|-----------------------------------------------------------------------|--------------------------------|--------------------------------------------------------|
| SRA_BAD_BIOSAMPLE_ID                    | BioSample accession appears to be invalid                             | no                             | Double check your input file  |
| SRA_FAIL_TO_DOWNLOAD_ALL                | ALL of a BioSample's run accessions fail prefetch and/or fasterq-dump | no (unless internet outage)    | If a few samples do this, it's probably because they are not actually paired-end Illumina FASTQs. If ALL samples fail like this, make sure you are able to access SRA using sra-tools: Try running the pipeline on just `SAMEA104362172` which is known to run to the end without issue. |
| SRA_FAIL_TO_DOWNLOAD_SOME               | ≥1 of a BioSample's run accessions fail prefetch and/or fasterq-dump  | yes                            | This can probably be safely ignored, and will not be an error unless `pull.fail_on_invalid` is true or if all run accessions fail. |
| SRA_ONE_FASTQ_ALL                       | ALL of a BioSample's run accessions have only one fastq               | no                             | Unless there is a way to split the FASTQ into two paired-end FASTQs, this sample cannot be used. |
| SRA_ONE_FASTQ_SOME                      | ≥1 of a BioSample's run accessions have only one fastq                | yes                            | This can probably be safely ignored, and will not be an error unless `pull.fail_on_invalid` is true or if all run accessions fail. |



<!--- 
| DECONTAMINATION_NOTHING_LEFT            | Comparing the number of reads in your FASTQ before and after decontamination indicates that the vast majority of it was contamination | yes, but the sample is suspect | Your sample was heavily contaminated! If your sample started out large enough, there might be enough data left to continue, which you can attempt with AAAAAAA.                           |
| EARLYQC_LOW_MEDIAN_COVERAGE             | TBProfiler detected your sample has a median coverage below AAAAAAAAAAAAAAAAAA                                                        | yes, but the sample is suspect | It's very likely that your sample would be filtered out by later coverage checks even if this check was skipped. If you wish to continue attempting to use it anyway, adjust AAAAAAAAAAA  | 
| TREE_TOO_MANY_LOW_COVERAGE_SITES        | 
  --->
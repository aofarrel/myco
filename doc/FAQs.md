# FAQs

## General
### Where do samples get dropped?
|   | pipeline                         | task                     | situation                                                                                                                       | can this filter be disabled?            | can be made a fatal error instead of a silent filter? |
|---|----------------------------------|--------------------------|---------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------|-------------------------------------------------------|
|   | myco_sra                         | get_sample_IDs           | BioSample accession appears to be invalid                                                                                       | no                                      | no                                                    |
|   | myco_sra                         | pull                     | ALL of a BioSample's run accessions fail prefetch and/or fasterq-dump                                                           | no                                      | yes, via `pull.fail_on_invalid`                       |
|   | myco_sra                         | pull                     | ≥1 of a BioSample's run accessions fail prefetch and/or fasterq-dump                                                            | yes, disabled by default                | yes, via `pull.fail_on_invalid`                       |
|   | myco_sra                         | pull                     | ALL of a BioSample's run accessions have only one fastq                                                                         | no                                      | yes, via `pull.fail_on_invalid`                       |
|   | myco_sra                         | pull                     | ≥1 of a BioSample's run accessions have only one fastq                                                                          | yes, disabled by default                | yes, via `pull.fail_on_invalid`                       |
|   | myco_sra, myco_raw               | decontam_each_sample     | sample takes ≥ `timeout_decontam_part1` minutes to map to the decontamination reference via `clockwork map_reads`               | yes, via `timeout_decontam_part1` = 0   | yes, via `variant_call_each_sample.crash_on_timeout`  |
|   | myco_sra, myco_raw               | decontam_each_sample     | sample takes ≥ `timeout_decontam_part2` minutes to map to the decontamination reference via `clockwork remove_contam`           | yes, via `timeout_decontam_part2` = 0   | yes, via `variant_call_each_sample.crash_on_timeout`  |
|   | myco_sra, myco_raw               | early_qc                 | the percentage (as float between 0 and 1) of reads in this sample's fqs below a quality score of 30 is > `early_qc_cutoff_q30`  | yes, disabled by default                | no                                                    |
|   | myco_sra, myco_raw, myco_cleaned | variant_call_each_sample | sample takes ≥ `timeout_variant_caller` minutes to map to the decontamination reference via `clockwork variant_call_one_sample` | yes, via `timeout_variant_caller` = 0   | yes, via `variant_call_each_sample.crash_on_timeout`  |
|   | myco_sra, myco_raw, myco_cleaned | variant_call_each_sample | non-timeout error in `clockwork variant_call_one_sample`                                                                        | no                                      | yes, via `variant_call_each_sample.crash_on_error`    |
|   | myco_sra, myco_raw, myco_cleaned | trees.cat_diff_files     | porportion of low coverage sites in a sample's diff file ≥ `diff_max_low_cov_pct_per_sample`                                             | yes, via `max_low_coverage_sites` = 1.0 | no                                                    |


Miscellanous notes:
* SRA data failing prefetch or fasterq-dump (from sra-tools) usually means that data is corrupt. If you find ENA data on SRA that looks it ought not be corrupt, but can't be downloaded, try ENABrowserTools or [my WDlization of it](https://github.com/aofarrel/enaBrowserTools-wdl)
* These timers apply to the same WDL task but are for different processes within that task -- `timeout_decontam_part1` is 20 and `timeout_decontam_part1` is 15, and a sample spends 19 minutes mapping plus another 14 minutes finishing the decontamination process, it will *not* be filtered out.
* clockwork's variant caller failing could mean it ran out of memory or was passed bad data, so the pipeline doesn't crash when this happens unless `variant_call_each_sample.crash_on_error` is set to true. Erroring out in the decontamination step is more indicative a serious issue, so a non-timeout error in that step will crash the pipeline.

### Where does data get filtered?
* If a FASTQ is above `subsample_cutoff` MB, it will get downsampled by seqtk
* When diff files are created, two types of regions will be masked:
  * regions that tend to be frequently masked when dealing with TB are masked (note: upstream outputs such as the VCF do not get masked)
  * regions which have called a variant but that variant has less coverage than 

### Why are there so many places where samples get dropped?
Myco was originally designed with the goal of analyzing as much TB SRA data as we could relatively quickly and cheaply. About 94% of MTBC BioSamples which are tagged as containing paired Illumina data pass myco_sra default's filters. An additional 3% of MTBC BioSamples are either completely invalid (fails fasterq-dump due to a database migration error that happended years ago, length of quality score doesn't match length of nucleotide strings, etc) or can't be used by clockwork (not actually paired Illumina reads, etc). The remainder of what gets filtered out seems to be extremely heavily contaminated and/or low overall coverage. As such, we're reasonably confident that myco_sra's default filters are a good tradeoff for the realities of SRA's data. Of course, you might want to cast a wider net than we did, or you might not be using SRA data at all, so all of these filters except for this-will-100%-break-the-pipeline-if-you-let-this-through ones can be configured.


### What if I want to use a different reference genome?
This isn't officially supported due to TBProfiler, UShER, and clockwork each needing specific reference genomes:
* TBProfiler's reference genome must be *exactly* the same as the one you called variants upon, if you're running TBProfiler on bams
* UShER has limitations on how long a chromosome's name can be
* clockwork's decontamination is designed with their specific decontamination reference in mind
That being said, old versions of myco used clockwork reference prepare to prepare the TB genome, and you could hack that passing-reference-genomes-around functionality to use your own custom reference genomes if you're confident. The latest version of myco_raw and myco_sra that used clockwork reference prepare was [4.1.3](https://github.com/aofarrel/myco/releases/tag/4.1.3), so that's a good place to start. The latest version of myco that required reference genomes for TBProfiler and/or UShER was [4.2.0](https://github.com/aofarrel/myco/releases/tag/4.2.0).


## Common warnings/errors
### trees/tree_nine/cat_diff_files fails with `Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE' but got: 'local-disk 0 SSD'`
This means you enabled the phylogenetic tree task, but none of your samples actually made it to the phylogenetic tree building task. See "Where do samples get dropped?" for more information. Due to the way Terra's UI works, it was decided it is easier to keep this as a visible error to make troubleshooting easier.
# Running myco

If you are not familiar with running WDLs, [please see this guide](https://github.com/ucsc-cgp/training-resources/blob/main/WDL/running_a_wdl.md).

## Inputs
| Name                        | Type                | Default | Info                                                                                                                                                                                                                                                                        |
|---------------------------  |-----------------    |-------- |---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  |
| **biosample_accessions**    | File                | n/a     | *(myco_sra.wdl only)* Text file listing BioSample accessions to pull reads from. Each line should have only one accession. SRS, SAM, ERS, and integer inputs are all supported. [Here's a bunch of sample files you can use!](https://github.com/aofarrel/SRANWRP/tree/main/inputs/quick_tests)   |
| **paired_fastqs**           | Array[Array[File]]  | n/a     | *(myco.wdl only)* Nested array of FASTQ files, where each internal array represents one sample's FASTQ files. Must be paired-end Illumina reads. [Here's a bunch of sample files you can use!](https://github.com/aofarrel/SRANWRP/tree/main/inputs/quick_tests) |
| less_scattering             | Bool                | false   | *(myco_sra.wdl only)* Set to `true` to prevent the decontamination process from executing as a scattered task. **It is recommended you leave this as false.**                   |
| min_coverage                | Int                 | 10      | Minimum coverage for a variant to be considered. Only affects the diff file, not the VCF.     |
| subsample_cutoff            | Int                 | -1      | If not set to -1, any FASTQ larger than this in value in MB will be downsampled with `seqtk`.     |
| subsample_seed              | Int                 | 1965    | Seed to use when downsampling with `seqtk`. Unused if `subsample_cutoff` is -1.   |
| tar_fqs                     | Bool                | false   | *(myco_sra.wdl only)* Set to `true` to tarball fastq files after downloading them from SRA. **Required if less_scattering = true**   |
| typical_tb_masked_regions   | File                | n/a     | BED file of ornery regions to mask. We recommend using [this mask file](https://github.com/iqbal-lab-org/cryptic_tb_callable_mask/blob/43ec21319209b23f648f32e4868bdf07cf09f2a0/R00000039_repregions.bed). Only affects the diff file, not the VCF.                      |

Note that *subsample_cutoff* and *subsample_seed* are workflow-level inputs in *myco.wdl* but are task-level inputs (of the pull task) in *myco_sra.wdl*. This is relevent only if you're using a JSON file to enter your inputs.

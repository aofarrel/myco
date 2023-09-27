# myco 🍄
myco is group of pipelines built for phylogenic analysis of the _Mycobacterium tuberculosis_ complex (MBTC). It builds upon existing tools such as [clockwork](https://github.com/iqbal-lab-org/clockwork) and [UShER](https://www.nature.com/articles/s41588-021-00862-7) to accomplish this task.

## Which workflow should I use?
Each version of myco largely only differs in how you are passing in FASTQ files. **In all cases, your FASTQs must be paired-end Illumina reads.**
* [myco_simple](https://qa.dockstore.org/workflows/github.com/aofarrel/myco/myco_simple) expects decontaminated and merged FASTQs
* [myco_raw](https://qa.dockstore.org/workflows/github.com/aofarrel/myco/myco_raw) expects FASTQs which have not been decontaminated and may or may not be merged
* [myco_sra](https://qa.dockstore.org/workflows/github.com/aofarrel/myco/myco_sra) expects a text file listing BioSample accessions you wish to pull FASTQs from

For more information please see [./docs/inputs.md](./doc/inputs.md).

## More information
* How to use WDL workflows: [UCSC's guide on running WDLs](https://github.com/ucsc-cgp/training-resources/blob/main/WDL/running_a_wdl.md)
* Full list of inputs: [inputs.md](./doc/inputs.md)
* Per-workflow readmes:
  * [myco_raw](./doc/myco_raw.md)
  * [myco_simple](./doc/myco_simple.md)
  * [myco_sra](./doc/myco_sra.md)
  * [wrapper_example](.doc/wrapper_example.md)

myco imports almost all of its code from other repos. Please see those specific repos for support with different parts of the myco pipeline:
* Downloading reads from SRA: [SRANWRP](https://github.com/aofarrel/SRANWRP)
* Decontamination and calling variants: [clockwork-wdl](https://github.com/aofarrel/clockwork-wdl)
* Turning VCFs into diff files: [parsevcf](https://github.com/lilymaryam/parsevcf)
* Building UShER, Taxonium, and Nextstrain/Auspice trees: [tree-nine](https://github.com/aofarrel/tree-nine)
* Determining why your sample spent 20 hours in the decontamination task: [FastQC-wdl](https://qa.dockstore.org/workflows/github.com/aofarrel/fastqc-wdl/fastqc)

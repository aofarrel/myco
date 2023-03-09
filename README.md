# myco 🍄
myco is a pipeline built for phylogenic analysis of the _Mycobacterium tuberculosis_ complex (MBTC). It builds upon existing tools such as [clockwork](https://github.com/iqbal-lab-org/clockwork) and [UShER](https://www.nature.com/articles/s41588-021-00862-7) to accomplish this task.

## Which workflow should I use?
Each version of myco is relatively similar. Where they differ is the sort of the inputs they expect. **In all cases, your FASTQs must be paired-end Illumina reads.**
* pairs of FASTQs which have been decontaminated and merged such that each sample has precisely two FASTQs associated with it: myco_cleaned
* pairs of FASTQs which have yet to be decontaminated or merged:
  * if each sample has its FASTQs in a single array: myco_raw
  * if each sample has its foward FASTQs in one array and reverse FASTQs in another array: [Decontam_And_Combine_One_Samples_Fastqs](https://dockstore.org/workflows/github.com/aofarrel/clockwork-wdl/Decontam_And_Combine_One_Samples_Fastqs), then myco_cleaned
* a list of SRA BioSamples whose FASTQs you'd like to use: myco_sra
* a list of SRA run accessions (ERR, SRR, DRR) whose FASTQs you'd like to use: [convert them to BioSamples](https://dockstore.org/workflows/github.com/aofarrel/SRANWRP/get_biosample_accessions_from_run_accessions:main?tab=info), then myco_sra

## More information
* Basic info on running myco: [running_myco.md](./doc/running_myco.md)
* Full list of inputs: [inputs.md](./doc/inputs.md)
* Overview of the pipeline: [pipline_overview.md](./doc/pipline_overview.md)

myco imports almost all of its code from other repos. Please see those specific repos for support with different parts of the myco pipeline:
* Downloading reads from SRA: [SRANWRP](https://github.com/aofarrel/SRANWRP)
* Decontamination and calling variants: [clockwork-wdl](https://github.com/aofarrel/clockwork-wdl)
* Turning VCFs into diff files: [parsevcf](https://github.com/lilymaryam/parsevcf)
* Building UShER and Taxonium trees: [usher-sampled-wdl](https://github.com/aofarrel/usher-sampled-wdl)


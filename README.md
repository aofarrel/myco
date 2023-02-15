# myco 🍄
myco is a pipeline built for phylogenic analysis of the _Mycobacterium tuberculosis_ complex (MBTC). It builds upon existing tools such as [clockwork](https://github.com/iqbal-lab-org/clockwork) and [UShER](https://www.nature.com/articles/s41588-021-00862-7) to accomplish this task.

myco imports almost all of its code from other repos. Please see those specific repos for support with different parts of the myco pipeline:
* Downloading reads from SRA: [SRANWRP](https://github.com/aofarrel/SRANWRP)
* Decontamination and calling variants: [clockwork-wdl](https://github.com/aofarrel/clockwork-wdl)
* Turning VCFs into diff files: [parsevcf](https://github.com/lilymaryam/parsevcf)
* Building UShER and Taxonium trees: [usher-sampled-wdl](https://github.com/aofarrel/usher-sampled-wdl)

For more information, please see [docs/readme.md](./docs/readme.md).

## Which workflow should I use?
If you already have a bunch of fastqs: myco.wdl  
If you want to pull fastqs from SRA: myco_sra.wdl 
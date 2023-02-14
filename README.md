# myco 🍄
myco is a pipeline built for phylogenic analysis of the _Mycobacterium tuberculosis_ complex (MBTC). It builds upon existing tools such as [clockwork](https://github.com/iqbal-lab-org/clockwork) and [UShER](https://www.nature.com/articles/s41588-021-00862-7) to accomplish this task.

myco imports almost all of its code from other repos. Please see those specific repos for support with different parts of the myco pipeline:
* Downloading reads from SRA: [SRANWRP](https://github.com/aofarrel/SRANWRP)
* Decontamination and calling variants: [clockwork-wdl](https://github.com/aofarrel/clockwork-wdl)
* Turning VCFs into diff files: [parsevcf](https://github.com/lilymaryam/parsevcf)
* Building UShER and Taxonium trees: [usher-sampled-wdl](https://github.com/aofarrel/usher-sampled-wdl)

## Which workflow should I use?
If you already have a bunch of fastqs: myco.wdl  
If you want to pull fastqs from SRA: myco_sra.wdl  

Note that pulling from SRA allows for more data validation than putting in fastqs directly. If using myco.wdl, please make sure your data:
* is Illumina paired-end data  
* is grouped per-sample  
* is actually tuberculosis  
* len(quality scores) = len(nucleotides) for every line  
* isn't huge — individual files over `subsample_cutoff` (default: 450 MB) will be downsampled, but keep an eye on the cumulative size of samples which have lots of small reads (ex: SAMEA968096)


## Running myco
The currently recommended way to run myco is on Terra, a cloud-computing resource which uses Google as its backend. However, it also runs as-is on local machines. See [running_myco.md](/doc/running_myco.md) for myco-specific instructions, and [this guide](https://github.com/ucsc-cgp/training-resources/blob/main/WDL/running_a_wdl.md) if you're new to running WDLs in general.
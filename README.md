# myco 🍄
myco is a pipeline built for phylogenic analysis of the _Mycobacterium tuberculosis_ complex (MBTC). It builds upon existing tools such as [clockwork](https://github.com/iqbal-lab-org/clockwork) and [UShER](https://www.nature.com/articles/s41588-021-00862-7) to accomplish this task.

myco imports almost all of its code from other repos. Please see those specific repos for support with different parts of the myco pipeline:
* Downloading reads from SRA: [SRANWRP](github.com/aofarrel/SRANWRP)
* Decontamination and calling variants: [clockwork-wdl](github.com/aofarrel/clockwork-wdl)
* Turning VCFs into diff files: [parsevcf](https://github.com/lilymaryam/parsevcf)


## Running myco
The currently recommended way to run myco is on Terra, a cloud-computing resource which uses Google as its backend. However, it also runs as-is on local machines. See [running_myco.md](./docs/running_myco.md) for instructions.
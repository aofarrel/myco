# Reference genome information
Reference genomes are handled automatically, but there are a few things the typical user should know.

As of release 4.3.0, all Docker images used by myco now come with their own copy of the reference genome by packaging the appropriate output of clockwork refprep. Refprep is my implementation of [clockwork's reference preparation standards](https://github.com/iqbal-lab-org/clockwork/wiki/Walkthrough-scripts-only#get-and-index-reference-genomes). For most tasks, this is just H37Rv, for the decontamination task, this also includes a decontamination reference.

**Because the decontamination reference is quite large (~12 GB), it will take some time to download from Docker Hub if it's not already present on your machine. This may cause the decontamination task to look "stuck" before it starts properly executing.**

## More on the decontamination reference
The decontamination reference I use is the same one clockwork v0.12.5 uses. This is *not* the same decontamination reference CDC uses, which you can use instead by flipping `decontam_use_CDC_varpipe_ref` to true, but I do not recommended doing so unless compliance with CDC's standards is a hard requirement

The following Docker images contain the two decontamination reference genomes I support for this pipeline.
* ashedpotatoes/clockwork-plus:v0.12.5.3-CRyPTIC
    * This is the default
    * Also contains clockwork-plus:v0.12.5
* ashedpotatoes/clockwork-plus:v0.12.5.3-CDC
    * Uses the same decontamination reference used by [CDC's NCHHSTP-DTBE-Varpipe-WGS workflow](https://github.com/CDCgov/NCHHSTP-DTBE-Varpipe-WGS)
    * Functional, but missing some contigs and metadata, so not recommended unless you must strictly replicate CDC results
    * Also contains clockwork-plus:v0.12.5
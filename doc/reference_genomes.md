# Reference genome information
Reference genomes are handled automatically, but there are a few things the typical user should know.

myco executes in Docker images, which are primarily used to bundle software, but can include any arbitrary files. As of release 4.3.0, all Docker images used by myco now come with their own copy of the H37Rv reference [NC_000962.3](https://www.ncbi.nlm.nih.gov/nuccore/448814763). The decontamination task `clckwrk_combonation.clean_and_decontam_and_check`, which may show on your logs as `fastp_decontam_check`, use a different Docker image that additionally contains a decontamination reference. Because the decontamination reference is quite large (~12 GB), the resulting Docker image is also quite large. **It may take several minutes to download this image Docker Hub if it's not already present on your machine. This may cause the decontamination task to look "stuck" before it starts properly executing.**

## More on the decontamination reference
The decontamination reference I use is the same one clockwork v0.12.5 uses. This is *not* the same decontamination reference CDC uses, which you can use instead by flipping `decontam_use_CDC_varpipe_ref` to true, but I do not recommended doing so unless compliance with CDC's standards is a hard requirement.

The following Docker images contain the two decontamination reference genomes I support for this pipeline.
* ashedpotatoes/clockwork-plus:v0.12.5.3-CRyPTIC
    * This is the default, as it's the reference specifically made for the clockwork's decontamination task, built by [clockwork's reference preparation workflow](https://github.com/iqbal-lab-org/clockwork/wiki/Walkthrough-scripts-only#get-and-index-reference-genomes)
    * [Contents can be viewed here](https://github.com/iqbal-lab-org/clockwork/blob/b3159f39cb417eb42711ade4a9cd08cd02dea7d2/scripts/download_tb_reference_files.pl)
    * Human decontamination reference based on CHM13
    * Docker image also contains clockwork v0.12.5 software
* ashedpotatoes/clockwork-plus:v0.12.5.3-CDC
    * Uses the same decontamination reference used by [CDC's NCHHSTP-DTBE-Varpipe-WGS workflow](https://github.com/CDCgov/NCHHSTP-DTBE-Varpipe-WGS)
    * Functional, but missing some contigs and metadata, so not recommended unless you must strictly replicate CDC results
    * Unknown provenance *(if you are familiar with this decontamination reference, please drop a PR with some information!)*
    * Human decontamination reference likely hg38
    * Docker image also contains clockwork v0.12.5 software
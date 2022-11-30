# Running myco

## Backends

### On Terra
These instructions assume you already have an account on Terra and [a billing project set up](https://support.terra.bio/hc/en-us/articles/360026182251-How-to-set-up-billing-in-Terra).
1. Go to myco's page on Dockstore
2. Click the button on the right hand side to import into Terra
3. Select which Terra workspace you wish to import into, or select a new one
4. Go to the workflow tab, select myco, and run your workflow

Note: It is also possible to run myco by simply copy-pasting the WDL via the Broad Methods Repository, but doing this will not transfer over git versioning, nor will your copy of the workflow keep up-to-date automatically.

### On a local machine (ie: your laptop)
[miniwdl](https://github.com/chanzuckerberg/miniwdl) generally works best for running myco locally, as it seems to handle local resources better than Cromwell (ie, doesn't crash the Docker daemon). Since this WDL uses mv commands on input files, you'll need to use the --copy-input-files option.

`miniwdl run myco.wdl --copy-input-files -i input/per_sample/quickmini.json`

Cromwell is what Terra uses on the backend, so if you want are debugging for Terra, I recommend using the Docktore CLI or latest version of Cromwell.

### On an HPC
Like most WDLs, myco uses a Docker image. However, many (not all) institutes do not allow Docker to run on their HPC systems for security reasons. Strictly speaking, [there are ways](https://docs.dockstore.org/en/stable/advanced-topics/docker-alternatives.html) to get around this limitation, but we cannot provide support for HPC users whose HPC administrators do not allow for running Docker images.

## Inputs
* samples: Text file listing BioSample accessions to pull reads from. Each line should have only one accession. SRS, SAM, and integer inputs are all supported.
* min_coverage: Minimum coverage required to show up in the final diff file
* typical_tb_masked_regions: Currently we're using [this mask file](https://github.com/iqbal-lab-org/cryptic_tb_callable_mask/blob/43ec21319209b23f648f32e4868bdf07cf09f2a0/R00000039_repregions.bed) but you could insert your own.
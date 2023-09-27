# myco_simple
myco_simple is the version of myco to use if you already have a bunch of fastqs, divided on a per-sample basis, and at least one of these is true:
* your fastqs are already decontaminated and cleaned
* you do not care about decontaminating/cleaning your fastqs
* you want to run a very basic version of myco as a test
* you want to compare the final output of a different version of myco versus what you would get without decontamination/fastp cleaning

## Notable inputs
All non-fastq inputs are documented here: [inputs.md](./inputs.md)

### FASTQs
You need your FASTQs as a nested array, where each inner array represents one sample. For example, let's say you the following samples in a google bucket located at gs://my-cool-bucket/fqs/
* SAMN02599053, consisting of SAMN02599053_SRR1173122_1.fq.gz, SAMN02599053_SRR1173122_2.fq.gz, SAMN02599053_SRR1173191_1.fq.gz, and SAMN02599053_SRR1173191_2.fq.gz
* SAMN13813990, consisting of SAMN13813990_SRR10869128_1.fq.gz and SAMN13813990_SRR10869128_2.fq.gz

You would input the following for **paired_fastq_sets**:

```
[["gs://my-cool-bucket/fqs/SAMN02599053_SRR1173122_1.fq.gz", "gs://my-cool-bucket/fqs/SAMN02599053_SRR1173122_2.fq.gz", "gs://my-cool-bucket/fqs/SAMN02599053_SRR1173191_1.fq.gz", "gs://my-cool-bucket/fqs/SAMN02599053_SRR1173191_2.fq.gz"], ["gs://my-cool-bucket/fqs/SAMN13813990_SRR10869128_1.fq.gz", gs://my-cool-bucket/fqs/SAMN13813990_SRR10869128_2.fq.gz"]]
```

The first array represents sample SAMN02599053. The second array represents sample SAMN13813990. You will note that SAMN02599053 has two pairs of fastqs, while SAMN13813990 only has one pair of fastqs -- this is fine!

## Full workflow process

### [1] clockwork variant_call_single
Based on clockwork variant_call_single, which itself combines samtools, cortex, and minos. For each sample, the output is a single VCF file and a BAM file.

### [2] (optional) Run FastQC on slow samples
If a sample times out in the variant calling step, it is usually due to an issue with the inputs. FastQC examines all inputs that timed out so you can see what might be going on.

### [3] VCF2diff -- Mask the outputs and optionally create diff files
When feeding outputs into UShER, we want to make use of diff files. But first, we perform a little bit of data processing -- it common for some regions of the TB genome to be masked. We want to avoid those problematic regions in our final output, as well as any regions without much coverage. This task cleans up our outputs and optionally creates a diff file, one per sample, which can be used to make some happy little trees.

### [4] (optional) Tree Nine -- Generate UShER, Taxonium, and NextStrain trees
If decorate_trees = true, and an input tree is passed in, each sample will be placed on the tree by UShER. The resulting tree will then be converted to Taxonium format, allowing it to be viewed in taxonium. NextStrain subtree JSONs will also be generated.
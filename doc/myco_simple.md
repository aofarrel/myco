# myco_simple
myco_simple is the version of myco to use if you already have a bunch of fastqs, divided on a per-sample basis. It assumes you do not want to decontaminate your fastqs at all -- perhaps they are non-tuberculosis mycobacteria or already decontaminated -- but it does offer the option of cleaning with fastp.

Each sample's read direction must be an individual gzipped file, e.g. SRR1173122_1.fq.gz and SRR1173122_2.fq.gz representing sample SRR1173122. Something like SRR1173122.gz which contains two fastqs would not work. Unlike more flexible versions of myco, the inputs to myco_simple **strictly** must be individually gzipped due to limitations of the variant caller.

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
When feeding outputs into UShER, we want to make use of diff files. But first, we perform a little bit of data processing -- it common for some regions of the TB genome to be masked. We want to avoid those problematic regions in our final output, as well as any regions without much coverage. This task cleans up our outputs and optionally creates a diff file, one per sample, which can be used to make some happy little phylogenetic trees using [Tree Nine](https://dockstore.org/workflows/github.com/aofarrel/tree_nine/tree_nine).

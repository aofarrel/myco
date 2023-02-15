# Running myco
myco is a WDL workflow. The currently recommended way to run myco is on Terra, a cloud-computing resource which uses Google as its backend. However, it also runs as-is on local machines.

If you are not familiar with running WDLs, [please see this guide](https://github.com/ucsc-cgp/training-resources/blob/main/WDL/running_a_wdl.md).

## How do I get my fastqs into myco?
First of all, please make sure your data:
* is Illumina paired-end data <sup>†</sup>   
* is grouped per-sample <sup>†</sup>     
* len(quality scores) = len(nucleotides) for every line <sup>†</sup>   
* is actually MTBC
* isn't huge — individual files over `subsample_cutoff` (default: 450 MB) will be downsampled, but keep an eye on the cumulative size of samples which have lots of small reads
  * it is okay to have more than two reads per sample -- where things get iffy is if you have 8 or more fastqs per sample (ex: SAMEA968096)

<sup>†</sup> myco_sra.wdl is able to detect these issues and will throw out those samples without erroring. myco.wdl is not able to detect these issues. Neither will check if your sample is actually MTBC or if a single sample has lots of small reads.

### Pulling reads from SRA
If your reads are on SRA, you should use myco_sra.wdl. Simply input the path to a text file that contains the BioSample accessions of the samples you are interested in (one BioSample accession per line). For example, if you have a file named my_accessions.txt at gs://my-cool-bucket, and my_accessions.txt looks like this:
```
SAMN02599053
SAMN13813990
```

Then you would input `gs://my-cool-bucket/my_accessions.txt` for **biosample_accessions** in myco_sra.wdl

### Direct input
If you already have a bunch of fastqs, you will want to input their URIs as a nested array, where each inner array represents one sample. For example, let's say you the following samples in a google bucket located at gs://my-cool-bucket/fqs/
* SAMN02599053, consisting of SAMN02599053_SRR1173122_1.fq.gz, SAMN02599053_SRR1173122_2.fq.gz, SAMN02599053_SRR1173191_1.fq.gz, and SAMN02599053_SRR1173191_2.fq.gz
* SAMN13813990, consisting of SAMN13813990_SRR10869128_1.fq.gz and SAMN13813990_SRR10869128_2.fq.gz

You would input the following for **paired_fastq_sets** in myco.wdl:

```
[["gs://my-cool-bucket/fqs/SAMN02599053_SRR1173122_1.fq.gz", "gs://my-cool-bucket/fqs/SAMN02599053_SRR1173122_2.fq.gz", "gs://my-cool-bucket/fqs/SAMN02599053_SRR1173191_1.fq.gz", "gs://my-cool-bucket/fqs/SAMN02599053_SRR1173191_2.fq.gz"], ["gs://my-cool-bucket/fqs/SAMN13813990_SRR10869128_1.fq.gz", gs://my-cool-bucket/fqs/SAMN13813990_SRR10869128_2.fq.gz"]]
```

The first array represents sample SAMN02599053. The second array represents sample SAMN13813990. You will note that SAMN02599053 has two pairs of fastqs, while SAMN13813990 only has one pair of fastqs -- this is fine!




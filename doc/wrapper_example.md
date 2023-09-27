# wrapper_example

wrapper_example.wdl, aka myco_simple_one_sample, is just a wrapper for myco_simple. It is designed for [Terra data tables](https://support.terra.bio/hc/en-us/articles/360025758392) specifically. A Terra data table being used as an input for wrapper_example might look like this:

| entity:sample_id | FASTQ_forward             | FASTQ_reverse             | decontaminated_fastq_1 | decontaminated_fastq_2 |
|------------------|---------------------------|---------------------------|------------------------|------------------------|
| sampleA          | raw_A_r1.fq               | raw_A_r2.fq               | decontam_A_r1.fq       | decontam_A_r2.fq       |
| sampleB          | raw_B1_r1.fq, rawB2_r1.fq | raw_B1_r2.fq, rawB2_r2.fq | decontam_B_r1.fq       | decontam_B_r2.fq       |
| sampleC          | raw_C_r1.fq               | raw_C_r2.fq               | decontam_C_r1.fq       | decontam_C_r2.fq       |
 
 
**This means that myco_simple_one_sample runs a single sample per workflow (ie, one-workflow-one-sample)** This is in contrast to the non-1samp version of myco_simple, which runs as a single workflow with individual tasks scattered upon multiple samples (ie, one-workflow-multiple-samples). 

Basically, if you want to run myco_simple on more than one sample at a time, you can either: 
* run multiple instances of this wrapper workflow (easiest if your data is in a Terra data table)
* run myco_simple (easiest if your data is in a JSON or you are copy-pasting your fastqs directly into Terra's UI)

Everything in [doc/myco_simple.md](./myco_simple.md) also applies to this wrapper, with the exception of the FASTQs part.
# myco 
myco is group of pipelines built for phylogenic analysis of the _Mycobacterium tuberculosis_ complex (MBTC). It builds upon existing tools such as [clockwork](https://github.com/iqbal-lab-org/clockwork) and [UShER](https://www.nature.com/articles/s41588-021-00862-7) to accomplish this task.

In an amusing repeat of [somewhat questionable naming decisions made in 1896](https://wwwnc.cdc.gov/eid/article/14/3/et-1403_article), myco should not be confused with the similiarly-named fungal pathogen pipeline [MycoSNP](https://github.com/CDCgov/mycosnp-nf).

## Which workflow should I use?
Each version of myco largely only differs in how you are passing in FASTQ files. **In all cases, your FASTQs must be paired-end Illumina reads.**
* [myco_simple](https://qa.dockstore.org/workflows/github.com/aofarrel/myco/myco_simple) expects decontaminated, gzipped FASTQs 
* [myco_raw](https://qa.dockstore.org/workflows/github.com/aofarrel/myco/myco_raw) expects FASTQs which have not been decontaminated and may or may not be merged
* [myco_sra](https://qa.dockstore.org/workflows/github.com/aofarrel/myco/myco_sra) expects a text file listing BioSample accessions you wish to pull FASTQs from

For more information please see [./docs/inputs.md](./doc/inputs.md).

## More information
* How to use WDL workflows: [UCSC's guide on running WDLs](https://github.com/ucsc-cgp/training-resources/blob/main/WDL/running_a_wdl.md)
* Full list of inputs: [inputs.md](./doc/inputs.md)
* Per-workflow readmes:
  * [myco_raw](./doc/myco_raw.md)
  * [myco_simple](./doc/myco_simple.md)
  * [myco_sra](./doc/myco_sra.md)
  * [wrapper_example](.doc/wrapper_example.md)

myco imports almost all of its code from other repos. Please see those specific repos for support with different parts of the myco pipeline:
* Downloading reads from SRA: [SRANWRP](https://github.com/aofarrel/SRANWRP)
* Decontamination and calling variants: [clockwork-wdl](https://github.com/aofarrel/clockwork-wdl)
* Turning VCFs into diff files: [parsevcf](https://github.com/lilymaryam/parsevcf)
* Building UShER, Taxonium, and Nextstrain/Auspice trees: [tree-nine](https://github.com/aofarrel/tree-nine)
* Determining why your sample spent 20 hours in the decontamination task: [FastQC-wdl](https://qa.dockstore.org/workflows/github.com/aofarrel/fastqc-wdl/fastqc)

## Citations

#### clockwork
> Hunt, Martin, Brice Letcher, Kerri M. Malone, Giang Nguyen, Michael B. Hall, Rachel M. Colquhoun, Leandro Lima, et al. “Minos: Variant Adjudication and Joint Genotyping of Cohorts of Bacterial Genomes.” Genome Biology 23, no. 1 (December 2022): 147. https://doi.org/10.1186/s13059-022-02714-x.

#### Cortex
> Iqbal, Zamin, Mario Caccamo, Isaac Turner, Paul Flicek, and Gil McVean. “De Novo Assembly and Genotyping of Variants Using Colored de Bruijn Graphs.” Nature Genetics 44, no. 2 (February 2012): 226–32. https://doi.org/10.1038/ng.1028.

#### fastp
> Chen, Shifu. “Ultrafast One‐pass FASTQ Data Preprocessing, Quality Control, and Deduplication Using Fastp.” iMeta 2, no. 2 (May 2023): e107. https://doi.org/10.1002/imt2.107.

#### goleft
> https://github.com/brentp/goleft

#### matUtils
> McBroome, Jakob, Bryan Thornlow, Angie S. Hinrichs, Alexander Kramer, Nicola De Maio, Nick Goldman, David Haussler, Russell Corbett-Detig, and Yatish Turakhia. “A Daily-Updated Database and Tools for Comprehensive Sars-Cov-2 Mutation-Annotated Trees.” Molecular Biology and Evolution 38, no. 12 (December 9, 2021): 5819–24. https://doi.org/10.1093/molbev/msab264.

#### minimap2
> Li, Heng. “Minimap2: Pairwise Alignment for Nucleotide Sequences.” Edited by Inanc Birol. Bioinformatics 34, no. 18 (September 15, 2018): 3094–3100. https://doi.org/10.1093/bioinformatics/bty191.

#### samtools
> Danecek, Petr, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, et al. “Twelve Years of SAMtools and BCFtools.” GigaScience 10, no. 2 (January 29, 2021): giab008. https://doi.org/10.1093/gigascience/giab008.

#### seqtk
> https://github.com/lh3/seqtk

#### TBProfiler (version 4.4.2, database 2023-Mar-27)
> Phelan, Jody E., Denise M. O’Sullivan, Diana Machado, Jorge Ramos, Yaa E. A. Oppong, Susana Campino, Justin O’Grady, et al. “Integrating Informatics Tools and Portable Sequencing Technology for Rapid Detection of Resistance to Anti-Tuberculous Drugs.” Genome Medicine 11, no. 1 (December 2019): 41. https://doi.org/10.1186/s13073-019-0650-x.

Note that some versions of this pipeline specifically uses Thiagen's fork of TBProfiler. The pinned version of this fork that I use itself uses TBProfiler version 4.4.2, database 2023-Jan-19.

> Libuit, Kevin G., Emma L. Doughty, James R. Otieno, Frank Ambrosio, Curtis J. Kapsak, Emily A. Smith, Sage M. Wright, et al. 2023. “Accelerating Bioinformatics Implementation in Public Health.” Microbial Genomics 9 (7). https://doi.org/10.1099/mgen.0.001051.

#### Trimmomatic
> Bolger, Anthony M., Marc Lohse, and Bjoern Usadel. “Trimmomatic: A Flexible Trimmer for Illumina Sequence Data.” Bioinformatics 30, no. 15 (August 1, 2014): 2114–20. https://doi.org/10.1093/bioinformatics/btu170.

#### UShER
> Turakhia, Yatish, Bryan Thornlow, Angie S. Hinrichs, Nicola De Maio, Landen Gozashti, Robert Lanfear, David Haussler, and Russell Corbett-Detig. “Ultrafast Sample Placement on Existing tRees (UShER) Enables Real-Time Phylogenetics for the SARS-CoV-2 Pandemic.” Nature Genetics 53, no. 6 (June 2021): 809–16. https://doi.org/10.1038/s41588-021-00862-7.

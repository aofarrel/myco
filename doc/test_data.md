# test data


## contamination
### contaminated in silico
#### TB + 10% Haemophilus influenzae
```
["gs://ucsc-pathogen-genomics-public/tb/fq/contaminated/insilico/HIB/ERR2510523_1.10pct.DRR235028_1.fq.gz",
"gs://ucsc-pathogen-genomics-public/tb/fq/contaminated/insilico/HIB/ERR2510523_2.10pct.DRR235028_2.fq.gz"]
```

### TB + 10% HIV
TODO

#### TB + 10% SARS-CoV-2
TODO


### contaminated in vivo(?)

### pure contaminants
```
["gs://ucsc-pathogen-genomics-public/NTM/M_avium_complex_sp/SAMN35100092_SRR24581793_1.fastq", 
"gs://ucsc-pathogen-genomics-public/NTM/M_avium_complex_sp/SAMN35100092_SRR24581793_2.fastq"],
["gs://ucsc-pathogen-genomics-public/NTM/M_kansasii/SAMN05182615_SRR3665538_1.fastq", 
"gs://ucsc-pathogen-genomics-public/NTM/M_kansasii/SAMN05182615_SRR3665538_2.fastq"],
["gs://ucsc-pathogen-genomics-public/NTM/M_kiyosense/SAMD00505764_DRR392335_1.fastq",
"gs://ucsc-pathogen-genomics-public/NTM/M_kiyosense/SAMD00505764_DRR392335_2.fastq"]
```

## "multistrain"
This is a set of samples that TBProfiler assigns multiple conflicting strains to.

## tree outliers
SAMEA838047
SAMEA4744414
SAMEA787734



### downsampled
	[ "gs://fc-1297615c-83b2-4119-a247-798ad975fbee/submissions/5fd1bf03-53ad-4395-b678-c798ab20c967/BIOSAMP_YOINK_NOFILE/03cc6292-f84e-41c9-a77a-7307a8a109e2/call-pull/shard-0/glob-db248e3bce81b54f5ef521878fe9e9de/SAMEA838047_ERR017800_1.fastq", "gs://fc-1297615c-83b2-4119-a247-798ad975fbee/submissions/5fd1bf03-53ad-4395-b678-c798ab20c967/BIOSAMP_YOINK_NOFILE/03cc6292-f84e-41c9-a77a-7307a8a109e2/call-pull/shard-0/glob-db248e3bce81b54f5ef521878fe9e9de/SAMEA838047_ERR017800_2.fastq", "gs://fc-1297615c-83b2-4119-a247-798ad975fbee/submissions/5fd1bf03-53ad-4395-b678-c798ab20c967/BIOSAMP_YOINK_NOFILE/03cc6292-f84e-41c9-a77a-7307a8a109e2/call-pull/shard-0/glob-db248e3bce81b54f5ef521878fe9e9de/SAMEA838047_ERR019874_1.fastq", "gs://fc-1297615c-83b2-4119-a247-798ad975fbee/submissions/5fd1bf03-53ad-4395-b678-c798ab20c967/BIOSAMP_YOINK_NOFILE/03cc6292-f84e-41c9-a77a-7307a8a109e2/call-pull/shard-0/glob-db248e3bce81b54f5ef521878fe9e9de/SAMEA838047_ERR019874_2.fastq", "gs://fc-1297615c-83b2-4119-a247-798ad975fbee/submissions/5fd1bf03-53ad-4395-b678-c798ab20c967/BIOSAMP_YOINK_NOFILE/03cc6292-f84e-41c9-a77a-7307a8a109e2/call-pull/shard-0/glob-db248e3bce81b54f5ef521878fe9e9de/SAMEA838047_ERR026644_1.fastq", "gs://fc-1297615c-83b2-4119-a247-798ad975fbee/submissions/5fd1bf03-53ad-4395-b678-c798ab20c967/BIOSAMP_YOINK_NOFILE/03cc6292-f84e-41c9-a77a-7307a8a109e2/call-pull/shard-0/glob-db248e3bce81b54f5ef521878fe9e9de/SAMEA838047_ERR026644_2.fastq" ]
	
	[ "gs://fc-1297615c-83b2-4119-a247-798ad975fbee/submissions/5fd1bf03-53ad-4395-b678-c798ab20c967/BIOSAMP_YOINK_NOFILE/03cc6292-f84e-41c9-a77a-7307a8a109e2/call-pull/shard-1/glob-db248e3bce81b54f5ef521878fe9e9de/SAMEA4744414_ERR2653228_1.fastq", "gs://fc-1297615c-83b2-4119-a247-798ad975fbee/submissions/5fd1bf03-53ad-4395-b678-c798ab20c967/BIOSAMP_YOINK_NOFILE/03cc6292-f84e-41c9-a77a-7307a8a109e2/call-pull/shard-1/glob-db248e3bce81b54f5ef521878fe9e9de/SAMEA4744414_ERR2653228_2.fastq" ]
	
	[ "gs://ucsc-pathogen-genomics-public/tb/fq/SAMEA787734_ERR040138_1.fastq", "gs://ucsc-pathogen-genomics-public/tb/fq/SAMEA787734_ERR040138_2.fastq" ]
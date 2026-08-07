# Test data you may be interested in (poorly organized notes)

You should probably just use whatever's in the checker data table.

SAMEA838047

## contamination
### contaminated in silico
#### TB + 10% Haemophilus influenzae
```
["gs://ucsc-pathogen-genomics-public/tb/fq/contaminated/insilico/HIB/ERR2510523_1.10pct.DRR235028_1.fq.gz",
"gs://ucsc-pathogen-genomics-public/tb/fq/contaminated/insilico/HIB/ERR2510523_2.10pct.DRR235028_2.fq.gz"]
```

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
This is a set of samples that TBProfiler assigns multiple conflicting strains to and are referenced in the Tree Nine repo

## tree outliers
SAMEA838047  
SAMEA4744414  
SAMEA787734  

### downsampled
	
	[ "gs://ucsc-pathogen-genomics-public/tb/fq/SAMEA787734_ERR040138_1.fastq", "gs://ucsc-pathogen-genomics-public/tb/fq/SAMEA787734_ERR040138_2.fastq" ]

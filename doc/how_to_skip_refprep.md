# How to skip refprep
Refprep is my implementation of [clockwork's reference preparation standards](https://github.com/iqbal-lab-org/clockwork/wiki/Walkthrough-scripts-only#get-and-index-reference-genomes). It is generally the first task any version of myco will be running and takes about an hour total to do the following:
1. Download TB reference files
2. Index the decontamination reference
3. Index the H37Rv reference

This is a deterministic<sup>†</sup> subworkflow, and several WDL executors (including Terra) allow for cacheing of previous workflow outputs, so this process usually only runs once. If you are using a backend/executor that doesn't support call cacheing, you can skip this process by letting refprep run once, then inputting the following:
* ClockworkRefPrepTB.bluepeter__tar_tb_ref_raw
* ClockworkRefPrepTB.bluepeter__tar_indexd_dcontm_ref
* ClockworkRefPrepTB.bluepeter__tar_indexd_H37Rv_ref

These files are too large for me to provide, but here's the structure of these tars for reference:

### ClockworkRefPrepTB.bluepeter__tar_tb_ref_raw
```
Ref.download.tar
    ├── NC_000962.1.fa
    ├── NC_000962.2.fa
    ├── NC_000962.3.fa
    ├── remove_contam.fa.gz
    └── remove_contam.tsv
```
### ClockworkRefPrepTB.bluepeter__tar_indexd_dcontm_ref
```
 Ref.remove_contam.tar
  ├── ref.fa
  ├── ref.fa.fai
  ├── ref.fa.minimap2_idx
  └── remove_contam_metadata.tsv
```
### ClockworkRefPrepTB.bluepeter__tar_indexd_H37Rv_ref
```
 Ref.H37Rv.tar
  ├── ref.fa
  ├── ref.fa.fai
  ├── ref.fa.minimap2_idx
  └── ref.k31.ctx
```

<sup>†</sup> If there is significant update to the TB reference that gets pulled by this script, that would change the output of refprep.

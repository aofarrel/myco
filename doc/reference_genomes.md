> [!IMPORTANT]  
> You are currently on an outdated branch of myco that exists solely for reproducing published results. It is HIGHLY recommended you use [a more recent version](https://github.com/aofarrel/myco) in order to take advantage of new updates to clockwork, TBProfiler, and other dependencies.

# Reference genome information
As of release 4.3.0, all Docker images used by myco now come with their own copy of the reference genome by packaging the appropriate output of clockwork refprep. Refprep is my implementation of [clockwork's reference preparation standards](https://github.com/iqbal-lab-org/clockwork/wiki/Walkthrough-scripts-only#get-and-index-reference-genomes). Prior to release 4.3.0, refprep had to run in each version of myco unless the user can use call cacheing or provided "bluepeter" inputs, which could be confusing for users, so nowadays I've packaged the outputs myself. Refprep takes about an hour total to do the following:
1. Download TB reference files
2. Index the decontamination reference
3. Index the H37Rv reference

This is a deterministic<sup>†</sup> subworkflow, and several WDL executors (including Terra) allow for cacheing of previous workflow outputs, so this process usually only ran once pre-4.3.0. If you are using a backend/executor that doesn't support call cacheing, you could skip this process by letting refprep run once, then inputting the following:
* ClockworkRefPrepTB.bluepeter__tar_tb_ref_raw
* ClockworkRefPrepTB.bluepeter__tar_indexd_dcontm_ref
* ClockworkRefPrepTB.bluepeter__tar_indexd_H37Rv_ref

These files are too large for me to provide on GitHub, but here's the structure of these tars for reference:

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

<sup>†</sup> If there is significant update to the TB reference that gets pulled by this script, that would change the output of refprep and the contents of the Docker images would need to be updated.

# Errors

## Downloading from SRA
### Error 1
If you set fail_on_invalid == True and fasterq-dump returns an odd number of files that is not 3, then the WDL task will return 1. When fasterq-dump returns 3 files, it is usually one file we can throw out plus two valid read files, so even if fail_on_invalid == True the code will not fail on the 3 file case.

### Error 3
The error comes from fasterq-dump itself, not the WDL. Error 3 has been observed with SRR1180764 and SRR1180610. These accessions are from the same mouse, and it appears 0% of their reads map to a tuberculosis reference genome. Most likely the error stems from this fact. Not sure why you can still get fastqs from SRA's website but not the CLI, but them's the breaks.
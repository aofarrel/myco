# myco
 Draft repository that will someday do things.

## Inputs
 Inputs are explained in the WDL's parameter_meta. Note that fastqs and samples are indexed by the order they are input -- put a given sample at the same index as its files, because the task that works on these things takes the dot-product of both arrays.

 ![diagram explaining that that string at index 0 of the sample array needs to correspond with the array of files at index 0 of the array of file arrays](/docs/fastqs_and_samples.png)

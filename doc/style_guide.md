# Style Guide (draft)

## Percentages/Ratios
User facing inputs and outputs that are a percentage are written as an integer or float between 0 and 100. Something that is 90% will be written as 90 or as 90.0, not 0.90
    ✔️ Int covstatsQC_max_percent_unmapped = 2 
    X Int covstatsQC_max_percent_unmapped = 0.2

## Units for disk size
Unless otherwise specified, user-facing storage units are in gigabytes. For example, variantcalling_addl_disk being 100 is interpreted as "add 100 GB to the disk calculation."


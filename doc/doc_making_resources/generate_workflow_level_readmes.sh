# WIP -- not currently used
#
# relies on https://github.com/Nicceboy/python-markdown-generator

python3 << CODE
from markdowngenerator import MarkdownGenerator

refprep_body = "Runs my implementation of [clockwork's reference preparation standards](https://github.com/iqbal-lab-org/clockwork/wiki/Walkthrough-scripts-only#get-and-index-reference-genomes). For more information on the specifics and how to skip this task if your backend doesn't support call cacheing, see [how_to_skip_refprep.md](./how_to_skip_refprep.md)"

decontam_body = "Based on [clockwork's decontamination process](https://github.com/iqbal-lab-org/clockwork/wiki/Walkthrough-scripts-only#decontaminate-the-reads), which runs clockwork map_reads and clockwork remove_contam in a single WDL task. The output is a group of decontaminated fastq files.

This step will also merge FASTQs if a single sample has more than one pair of FASTQs. For example, SAMN02599053 has four fastqs associated with it: 
* SAMN02599053_SRR1173122_1.fq.gz
* SAMN02599053_SRR1173122_2.fq.gz
* SAMN02599053_SRR1173191_1.fq.gz
* SAMN02599053_SRR1173191_2.fq.gz

The decontamination step will output a single pair: SAMN02599053_1.fastq and SAMN02599053_2.fastq"

varcall_body = "Based on clockwork variant_call_single, which itself combines samtools, cortex, and minos. For each sample, the output is a single VCF file and a BAM file."

fastqc_body = "If a sample times out in the decontamination or variant calling steps, it is usually due to an issue with the inputs. FastQC examines all inputs that timed out so you can see what might be going on."

diff_body = "When feeding outputs into UShER, we want to make use of diff files. But first, we perform a little bit of data processing -- it common for some regions of the TB genome to be masked. We want to avoid those problematic regions in our final output, as well as any regions without much coverage. This task cleans up our outputs and optionally creates a diff file, one per sample, which can be used to make some happy little trees."

tree_body = "If decorate_trees = true, and an input tree is passed in, each sample will be placed on the tree by UShER. The resulting tree will then be converted to Taxonium format, allowing it to be viewed in taxonium. NextStrain subtree JSONs will also be generated."


with open("../myco_raw.md", "w") as myco_raw:
	for line in lines:
		myco_raw.write(line)


with open("../myco_sra.md", "w") as myco_sra:
    pass

with open("../myco_cleaned.md", "w") as myco_cleaned:
    pass

CODE

echo "done"
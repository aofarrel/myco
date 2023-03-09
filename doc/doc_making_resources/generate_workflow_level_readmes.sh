# This is very much a Works On My Machine sort of script!
# This relies on https://github.com/Nicceboy/python-markdown-generator

python3 << CODE
from markdowngenerator import MarkdownGenerator

with open("../myco_raw.md", "w") as myco_raw:
	for line in lines:
		myco_raw.write(line)


with open("../myco_sra.md", "w") as myco_sra:
    pass

with open("../myco_cleaned.md", "w") as myco_cleaned:
    pass

CODE

echo "done"
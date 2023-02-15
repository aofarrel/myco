# This is very much a Works On My Machine sort of script!
# To run this elsewhere, point to your own version of womtool,
# change ghead to head (unless you're using Mac w/ coreutils),
# and pip3 install git+https://github.com/Nicceboy/python-markdown-generator

echo "grabbing inputs from myco_sra..."
java -jar /Applications/womtool-76.jar inputs myco_sra.wdl > raw.txt
echo "grabbing inputs from myco..."
java -jar /Applications/womtool-76.jar inputs myco.wdl >> raw.txt
echo "processing..."
sort raw.txt > sorted.txt
uniq sorted.txt > unique.txt
ghead -n -2 unique.txt | tail -n +2 | cut -c 9- > cleaned.txt

# I know this can probably be done with sed, but I don't like sed
python3 << CODE
import re
from markdowngenerator import MarkdownGenerator

def extract_wf_info(input_list):
	variables = []
	for line in input_list:
		this_variable = {"name": str(re.search('([a-z, A-Z, _, 1-9])+"\:', line).group(0)[:-2]),
						"type": str(re.search('\: \"([a-z, A-Z, _, 1-9])+\??', line).group(0)[3:]),
						"default": str(re.search("default = ([a-z, A-Z, 0-9])+", line).group(0)) 
							if bool(re.search("default = ([a-z, A-Z, 0-9])+", line)) 
							else ""}
		variables.append(this_variable)
	return variables

def extract_task_info(input_list):
	variables = []
	for line in input_list:
		this_variable = {"task": str(re.search('([a-z, A-Z, _, 1-9])+\.', line).group(0)[:-1]),
						"name": str(re.search('([a-z, A-Z, _, 1-9])+"\:', line).group(0)[:-2]),
						"type": str(re.search('\: \"([a-z, A-Z, _, 1-9])+\??', line).group(0)[3:]),
						"default": str(re.search("default = ([a-z, A-Z, 0-9])+", line).group(0)) 
							if bool(re.search("default = ([a-z, A-Z, 0-9])+", line)) 
							else ""}
		variables.append(this_variable)
	return variables

with open("cleaned.txt", "r") as f:
	task_level = []
	workflow_level = []
	for line in f.readlines():
		if bool(re.search("([a-z, A-Z, _])+\.", line)):
			task_level.append(line)
		else:
			workflow_level.append(line)

workflow_level = extract_wf_info(workflow_level)
task_level = extract_task_info(task_level)
runtime = []
not_runtime = []
for variable in task_level:
	if variable["name"] in ["addldisk", "cpu", "memory", "preempt"]:
		runtime.append(variable)
	else:
		not_runtime.append(variable)


with MarkdownGenerator(filename="inputs.md", enable_write=False) as doc:
	doc.addHeader(2, "Workflow-level inputs")
	doc.addTable(dictionary_list=workflow_level)
	doc.addHeader(2, "Task-level inputs")
	doc.addHeader(3, "Software settings")
	doc.addTable(dictionary_list=not_runtime)
	doc.addHeader(3, "Hardware settings")
	doc.addTable(dictionary_list=runtime)


CODE

rm sorted.txt unique.txt cleaned.txt
echo "done"
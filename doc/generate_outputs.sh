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

def strip_junk(string):
	return string.replace("&gt;", ">").replace("&quot;", "'").replace(": ", "")

def get_task(line):
	return str(re.search('([a-z, A-Z, _, 1-9])+\.', line).group(0)[:-1])

def get_name(line):
	return str(re.search('([a-z, A-Z, _, 1-9])+"\:', line).group(0)[:-2])

def get_type(line):
	return str(re.search('\: \"([a-z, A-Z, _, 1-9])+\??', line).group(0)[3:])

def get_default(line):
	if bool(re.search('default = .+?(?=\))', line)):
		return str(re.search('= .+?(?=\))', line).group(0).replace("(", "")).replace("= ", "")
	else:
		 return ""

def extract_wf_info(input_list):
	variables = []
	for line in input_list:
		this_variable = {"name": get_name(line),
						"type": get_type(line),
						"default": get_default(line)}
		variables.append(this_variable)
	return variables

def extract_task_info(input_list):
	variables = []
	for line in input_list:
		this_variable = {"task": get_task(line),
						"name": get_name(line),
						"type": get_type(line),
						"default": get_default(line)}
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
for input_variable in task_level:
	# runtime variables
	if input_variable["name"] == "addldisk":
		input_variable["description"] = "Additional disk size, in GB, on top of auto-scaling disk size."
		runtime.append(input_variable)
	elif input_variable["name"] == "disk_size":
		input_variable["description"] = "Disk size, in GB. Note that since cannot auto-scale as it cannot anticipate the size of reads from SRA."
		runtime.append(input_variable)
	elif input_variable["name"] == "cpu":
		input_variable["description"] = "Number of CPUs (cores) to request from GCP."
		runtime.append(input_variable)
	elif input_variable["name"] == "memory":
		input_variable["description"] = "Amount of memory, in GB, to request from GCP."
		runtime.append(input_variable)
	elif input_variable["name"] == "preempt":
		input_variable["description"] = "How many times should this task be attempted on a preemptible instance before running on a non-preemptible instance?"
		runtime.append(input_variable)
	else:
		not_runtime.append(input_variable)

# extract parameter_meta for workflow-level variables
parameter_meta = []
in_parameter_meta = False
with open("myco.wdl", "r") as myco:
	for line in myco:
		if line.startswith("\tparameter_meta"):
			in_parameter_meta = True
			continue
		elif line.startswith("\t}") and in_parameter_meta:
			break
		elif in_parameter_meta and not line.startswith("}"):
			this_parameter = {"name": re.search("\S.+?(?=\:)", line).group(0),
							"description": strip_junk(re.search('(?=\:).+', line).group(0).replace('\"', ""))}
			parameter_meta.append(this_parameter)
		else:
			continue

for input_variable in workflow_level:
	value = input_variable["name"]
	if value in ["biosample_accessions", "paired_fastq_sets"]:
		input_variable["description"] = "fastq input -- please see running_myco.md for more information"
	else:
		parameter = next((parameter for parameter in parameter_meta if parameter["name"] == value), None)
		input_variable["description"] = parameter["description"]

with MarkdownGenerator(filename="doc/inputs.md", enable_write=False) as doc:
	doc.writeTextLine("See /inputs/example_inputs.json for examples.")
	doc.addHeader(2, "Workflow-level inputs")
	doc.addTable(dictionary_list=workflow_level)
	doc.addHeader(2, "Task-level inputs")
	doc.addHeader(3, "Software settings")
	doc.writeTextLine("If you are on a backend that does not support call cacheing, you can use the 'bluepeter' inputs to skip the download of the reference genome.")
	doc.addTable(dictionary_list=not_runtime)
	doc.addHeader(3, "Hardware settings")
	doc.writeTextLine("A note on disk size: On GCP backends, disk size is treated as a maximum. If your task goes above that limit, it will fail.")
	doc.addTable(dictionary_list=runtime)

with open("doc/inputs.md", "r") as markdown:
	stripped = []
	for line in markdown:
		stripped.append(strip_junk(line))
with open("doc/inputs.md", "w") as markdown:
	for line in stripped:
		markdown.write(line)



CODE

rm raw.txt sorted.txt unique.txt cleaned.txt
echo "done"
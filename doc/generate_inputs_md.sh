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

filename_vars = [
			"out", 
			"contam_out_1", "contam_out_2", "counts_out",
			"done_file", "no_match_out_1", "no_match_out_2", 
			"outfile"
			]

def strip_junk(string):
	return string.replace("&gt;", ">").replace("&quot;", "'").replace(": ", "").replace("&#x27;", "\`")

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
	elif input_variable["name"] == "retries":
		input_variable["description"] = "How many times should we retry this task if it fails after it exhausts all uses of preemptibles?"
		runtime.append(input_variable)
	elif input_variable["name"] == "ssd":
		input_variable["description"] = "If true, use SSDs for this task instead of HDDs"
		runtime.append(input_variable)
	else:
		# not a runtime variable
		if input_variable["name"] in filename_vars:
			input_variable["description"] = "Override default output file name with this string"
			not_runtime.append(input_variable)
		elif input_variable["name"] == "histograms":
			input_variable["description"] = "Should coverage histograms be output?"
			not_runtime.append(input_variable)
		elif input_variable["name"] == "crash_on_timeout":
			input_variable["description"] = "If this task times out, should it stop the whole pipeline (true), or should we just discard this sample and move on (false)?"
			not_runtime.append(input_variable)
		elif input_variable["name"] == "subsample_cutoff":
			input_variable["description"] = "If a fastq file is larger than than size in MB, subsample it with seqtk (set to -1 to disable)"
			not_runtime.append(input_variable)
		elif input_variable["name"] == "subsample_seed":
			input_variable["description"] = "Seed used for subsampling with seqtk"
			not_runtime.append(input_variable)
		# decontamination 
		elif input_variable["name"] == "threads":
			input_variable["description"] = "Try to use this many threads for decontamination. Note that actual number of threads also relies on your hardware."
			not_runtime.append(input_variable)
		# variant caller
		elif input_variable["name"] == "debug":
			input_variable["description"] = "Do not clean up any files and be verbose"
			not_runtime.append(input_variable)
		elif input_variable["name"] == "mem_height":
			input_variable["description"] = "cortex mem_height option. Must match what was used when reference_prepare was run (in other words do not set this variable unless you are also adjusting the reference preparation task)"
			not_runtime.append(input_variable)
		else:
			input_variable["description"] = "" # need this or else the table is missing a column
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
	doc.addHeader(3, "Runtime attributes")
	doc.writeTextLine("These variables adjust runtime attributes, which includes hardware settings. See https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/ for more information.")
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
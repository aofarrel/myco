echo "grabbing inputs from myco_sra..."
#java -jar /Applications/womtool-76.jar inputs myco_sra.wdl > raw.txt
echo "grabbing inputs from myco..."
#java -jar /Applications/womtool-76.jar inputs myco.wdl >> raw.txt
echo "processing..."
sort raw.txt > sorted.txt
uniq sorted.txt > unique.txt
ghead -n -2 unique.txt | tail -n +2 | cut -c 9- > cleaned.txt # some shells can get away with head
# sed hurts my head, so we're going to cheat and use Python.
python3 << CODE
import re

def extract_info(input_list):
	vars = []
	types = []
	defaults = []
	for line in input_list:
		vars.append(re.search('([a-z, A-Z, _, 1-9])+"\:', line).group(0)[:-2])
		types.append(re.search("([A-Z])\w+\??", line).group(0))
		if bool(re.search("default = ([a-z, A-Z, 0-9])+", line)):
			defaults.append(re.search("default = ([a-z, A-Z, 0-9])+", line).group(0))
		else:
			defaults.append("")
	return [vars, types, defaults]

with open("cleaned.txt", "r") as f:
	task_level = []
	workflow_level = []
	for line in f.readlines():
		if bool(re.search("([a-z, A-Z, _])+\.", line)):
			task_level.append(line)
		else:
			workflow_level.append(line)

workflow_level_lists = extract_info(workflow_level)
i = 0
while i<len(workflow_level_lists[0]):
	print(workflow_level_lists[0][i])
	print(workflow_level_lists[1][i])
	print(workflow_level_lists[2][i])
	i += 1

# extract task name for the task_level before calling extract_info()


CODE


#$(pwd)/doc/inputs.md

#([A-Z])\w+\?? # grab types
#default = ([a-z, A-Z, 0-9])+ # grab defaults


# rm raw.txt sorted.txt unique.txt cleaned.txt
echo "done"
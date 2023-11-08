import sys

total_bp = 0
nr_added_nodes = 0
nr_total_nodes = 0
nr_total_links = 0

for line in sys.stdin:
	fields = line.strip().split()
	if fields[0] == "S":
		total_bp += len(fields[2])
		if ("sniffles" in line) or ("delly" in line) or ("svarp" in line):
			nr_added_nodes += 1
		nr_total_nodes += 1
	if fields[0] == "L":
		nr_total_links += 1

print('total sequence represented in graph [bp]: ' + str(total_bp))
print('number of added nodes (sniffles, delly, SVarp): ' + str(nr_added_nodes))
print('number of total nodes: ' + str(nr_total_nodes))
print('number of total links: ' + str(nr_total_links))

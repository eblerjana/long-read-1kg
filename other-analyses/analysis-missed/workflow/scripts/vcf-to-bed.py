import sys

for line in sys.stdin:
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	start = int(fields[1]) - 1
	end = start + len(fields[3])
	print('\t'.join([fields[0], str(start), str(end)]))

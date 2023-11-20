import sys

name = sys.argv[1]

for line in sys.stdin:
	if line.startswith('##'):
		print(line.strip())
		continue
	if line.startswith('#'):
		fields = line.strip().split()
		fields[-1] = name
		print('\t'.join(fields))
		continue
	print(line.strip())

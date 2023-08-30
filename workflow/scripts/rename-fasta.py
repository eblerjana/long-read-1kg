import sys

postfix = sys.argv[1]

for line in sys.stdin:
	if line.startswith('>'):
		fields = line.strip().split()
		fields[0] = fields[0] + '_' + postfix
		print(fields[0])
		continue
	print(line.strip())
		

import sys

n_bases = 0

for line in sys.stdin:
	if line.startswith(">"):
		continue
	for c in line.strip():
		if not c in ["n", "N", "\n"]:
			sys.stdout.write(c)
			n_bases += 1
print(n_bases)

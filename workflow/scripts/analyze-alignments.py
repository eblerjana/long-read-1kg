import sys
import re


def parse_cigar(cigar):
	parsed_cigar = re.findall(r'(\d+)([A-Z,=]{1})', cigar)
	return parsed_cigar


n_matches = 0
n_insertions = 0
n_deletions = 0
n_soft_clipped = 0
n_hard_clipped = 0


n_reads_no_gaps = 0
n_reads_only_small_gaps = 0

processed = 0

for line in sys.stdin:
	fields = line.strip().split()
	if processed % 100000 == 0:
		print('Processed ' + str(processed) + ' alignments.')
	for f in fields:
		if f.startswith("cg:Z:"):
			cigar = parse_cigar(f[5:])
			no_gaps = True
			only_small_gaps = True
			for c in cigar:
				length = int(c[0])
				if c[1] == "=":
					n_matches += length
				elif c[1] == "I":
					no_gaps = False
					if length > 50:
						only_small_gaps = False
					n_insertions += length
				elif c[1] == "D":
					no_gaps = False
					if length > 50:
						only_small_gaps = False
					n_deletions += length
				elif c[1] == "S":
					n_soft_clipped += length
				elif c[1] == "H":
					n_hard_clipped += length
			if no_gaps:
				n_reads_no_gaps += 1
			if only_small_gaps:
				n_reads_only_small_gaps += 1
	processed += 1


print('Matches: ' + str(n_matches))
print('Insertions: ' + str(n_insertions))
print('Deletions: ' + str(n_deletions))
print('Soft-clipped: ' + str(n_soft_clipped))
print('Hard-clipped: ' + str(n_hard_clipped))	
print('')
print('#reads with no gaps: ' + str(n_reads_no_gaps))
print('#reads with no gaps > 50bp ' + str(n_reads_no_gaps))	

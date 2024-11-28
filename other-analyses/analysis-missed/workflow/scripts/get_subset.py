import sys
import argparse

parser = argparse.ArgumentParser(prog='get_subset.py', description=__doc__)
parser.add_argument('--present', default = [], nargs='+', help="names of the columns that must be set to True.")
parser.add_argument('--absent', default = [], nargs='+', help="names of the columns that must be set to False.")
args = parser.parse_args()

present = args.present
absent = args.absent

present_ids = []
absent_ids = []

for line in sys.stdin:
	fields = line.strip().split()
	if line.startswith('ID'):
		for p in present:
			present_ids.append(fields.index(p))
		for a in absent:
			absent_ids.append(fields.index(a))

		continue

	# only keep line if all present ids are true
	# and all absent ids are false
	print_line = True
	for p in present_ids:
		if fields[p] == "False":
			print_line = False
			break
	for a in absent_ids:
		if fields[a] == "True":
			print_line = False
			break

	if print_line:
		print('\t'.join([fields[1], fields[2], fields[3], fields[0]]))	

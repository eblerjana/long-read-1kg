import sys
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def parse_intersection_file(filename, svtype):
	total_a = set([])
	total_b = set([])
	match_a = set([])
	match_b = set([])
	unique_a = set([])
	unique_b = set([])

	sample = filename.split('_')[-1].split('.')[0]

	callset_a = None
	callset_b = None

	for line in open(filename, 'r'):
		fields = line.strip().split()
		if line.startswith('ID'):
			callset_a = fields[9].split('_')[-1]
			callset_b = fields[10].split('_')[-1]
			continue

		if (fields[5] != svtype) and (svtype != 'ALL'):
			continue

		ids_a = set([a for a in fields[11].split(';') if a != 'nan'])
		ids_b = set([a for a in fields[12].split(';') if a != 'nan'])
		matched = ids_a and ids_b

		for i in ids_a:
			total_a.add(i)
			if matched:
				match_a.add(i)
			else:
				unique_a.add(i)

		for i in ids_b:
			total_b.add(i)
			if matched:
				match_b.add(i)
			else:
				unique_b.add(i)
	
	assert len(total_a) == len(match_a) + len(unique_a)
	assert len(total_b) == len(match_b) + len(unique_b)
	# comparison_SAGA-raw/results/SAGA-aug_PAV/intersection/intersection_SAGA-aug_PAV_HG00733.tsv		
	# chr1-1577126-DEL->s66173>s66195-85      chr1    1577122 1577208 86      DEL     60      nan     nan     True    True    chr1-1577126-DEL->s66173>s66195-85      chr1-1577123-DEL-85
	print('\t'.join([svtype, sample, callset_a, str(len(total_a)), str(len(match_a)), str(len(unique_a))]))
	print('\t'.join([svtype, sample, callset_b, str(len(total_b)), str(len(match_b)), str(len(unique_b))]))

	return sample, callset_a, len(match_a), len(unique_a), callset_b, len(match_b), len(unique_b)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='plot-intersections.py', description=__doc__)
	parser.add_argument('-o', '--output', required=True, help='Prefix of output files.')
	args = parser.parse_args()

	print('\t'.join(["vartype", "sample", "dataset", "all", "shared", "unique"]))
	filenames = [f for f in sys.stdin]

	with PdfPages(args.output) as pdf:
		for vartype in ['ALL', 'INS', 'DEL', 'OTHER', "COMPLEX", "INV", "DUP", "INVDUP", "CNV", "DUP:TANDEM"]:
			to_plot = {}
			samples = set([])
			c_a = None
			c_b = None
			for filename in filenames:
				sample, callset_a, match_a, unique_a, callset_b, match_b, unique_b = parse_intersection_file(filename.strip(), vartype)
				to_plot[(sample, callset_a)] = [match_a, unique_a]
				to_plot[(sample, callset_b)] = [match_b, unique_b]
				samples.add(sample)
				c_a = callset_a
				c_b = callset_b
			samples = sorted(list(samples))
			for c in [c_a, c_b]:
				fig = plt.figure()
				plt.title(c + ' ' + vartype)
				y_1 = [to_plot[(s, c)][0] for s in samples]
				y_2 = [to_plot[(s, c)][1] for s in samples]
				plt.xticks(rotation='vertical')
				plt.bar(samples, y_2, label = 'unique')
				plt.bar(samples, y_1, bottom = y_2, label = 'shared')
				plt.ylabel('number of variants')
				plt.legend()
				fig.tight_layout()
				pdf.savefig()
				plt.close()
			

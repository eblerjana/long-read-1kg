import sys
import argparse
import numpy as np
from matplotlib import pyplot as plt

def parse_tsv(filename):
	data = []
	for line in open(filename, 'r'):
		fields = line.split()
		# length of the longest branch of a bubble
		length = fields[7]
		data.append(int(length))
	return data


def parse_vcf(filename):
	data = []
	for line in open(filename, 'r'):
		if line.startswith('#'):
			continue
		fields = line.split()
		alleles = [len(fields[3])] + [len(a) for a in fields[4].split(',')]
		data.append(max(alleles))
	return data


def plot_histogram(datasets, labels, filename):
#	bins = np.histogram(np.hstack((d for d in datasets)), bins=100)[1]
	max_val = 0
	for data in datasets:
		for d in data:
			if d > max_val:
				max_val = d 
	for data, label in zip(datasets, labels):
		print(label)
		print('min:' + str(min(data)))
		print('max: ' + str(max(data)))
		plt.hist(data, bins=np.logspace(np.log10(10),np.log10(max_val + 10000), 200), alpha=0.4, label=label)
	plt.yscale('log')
	plt.xscale('log')
	plt.ylabel('Count')
	plt.xlabel('Bubble length (= length of longest path in bubble)')
	plt.legend(loc='upper right')
	plt.savefig(filename)

def run_plotting(filenames, labels, outname):
	datasets = []
	for f in filenames:
		if f.endswith('.tsv'):
			datasets.append(parse_tsv(f))
		elif f.endswith('.vcf'):
			datasets.append(parse_vcf(f))
		else:
			raise RuntimeError('Unknown file type of file ' + f)
	
	plot_histogram(datasets, labels, outname)
			


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='plot-histogram.py', description=__doc__)
	parser.add_argument('-files', metavar='FILE', nargs='+', required=True, help='File in TSV or VCF format.')
	parser.add_argument('-labels', metavar='LABELS', nargs='+', required=True, help='Labels. One per input file.')
	parser.add_argument('-outname', metavar='OUTNAME', default = 'histograms.pdf', help='name of the output file.')
	args = parser.parse_args()
	
	filenames = [f for f in args.files]
	labels = [l for l in args.labels]
	
	print(filenames, labels)
	if len(filenames) != len(labels):
		raise RuntimeError('Number of labels does not match number of files given.')
	
	run_plotting(filenames, labels, args.outname)

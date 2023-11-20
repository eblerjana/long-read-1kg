import sys
import argparse
import numpy as np
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(prog='mendelian-consistency.py', description=__doc__)
parser.add_argument('--min-count', metavar='MINCOUNT', type=int, required=True, help='Minimum number of carrier samples.')
parser.add_argument('--max-count', metavar='MAXCOUNT', type=int, required=True, help='Maximum number of carrier samples.')
parser.add_argument('--outname', metavar='OUTNAME', required=True, help='Name of the output PDF.')
args = parser.parse_args()

histogram = []
only_one_sample = 0
only_two_samples = 0
up_to_ten_samples = 0
more_than_ten_samples = 0
total_alleles = 0
lines_written = 0
max_val = 0

for line in sys.stdin:
	if line.startswith('#'):
		print(line.strip())
		continue
	fields = line.strip().split()
	nr_carriers = 0
	for genotype in fields[9:]:
		if '1' in genotype:
			nr_carriers += 1
	
	if nr_carriers == 1:
		only_one_sample += 1
	if nr_carriers == 2:
		only_two_samples += 1
	if nr_carriers <= 10:
		up_to_ten_samples += 1
	if nr_carriers > 10:
		more_than_ten_samples += 1

	if (nr_carriers >= args.min_count) and (nr_carriers <= args.max_count):
		print(line.strip()) 
		lines_written += 1
	total_alleles += 1
	histogram.append(nr_carriers)
	if nr_carriers > max_val:
		max_val = nr_carriers

sys.stderr.write('Number of alleles present in a single sample:\t' + str(only_one_sample) + '\n')
sys.stderr.write('Number of alleles present in two samples:\t' + str(only_two_samples) + '\n')
sys.stderr.write('Number of alleles present in up to 10 samples:\t' + str(up_to_ten_samples) + '\n')
sys.stderr.write('Number of alleles present in more than 10 samples:\t' + str(more_than_ten_samples) + '\n')
sys.stderr.write('Wrote ' + str(lines_written) + ' of ' + str(total_alleles) + ' to output VCF.' + '\n')

# plot histogram
plt.hist(histogram, bins=np.logspace(np.log10(1),np.log10(max_val + 10000), 200))
plt.yscale('log')
plt.xscale('log')
plt.ylabel('Count')
plt.xlabel('Number of carriers')
plt.savefig(args.outname)

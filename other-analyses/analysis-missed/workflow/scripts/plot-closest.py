import sys
import argparse
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt

def parse_tsv():
	data = []
	for line in sys.stdin:
		fields = line.strip().split()
		# distance to the closest feature in other BED
		distance = int(fields[-1])
		data.append(distance)
	return data



def plot_histogram(dataset, filename):

	with PdfPages(filename) as pdf:
		plt.figure()
		print('number of matches (distance=0): ' + str(dataset.count(0)))
		print('min distance:' + str(min(data)))
		print('max distance: ' + str(max(data)))
		plt.hist(dataset, bins=200, alpha=0.7)
		plt.yscale('log')
		plt.ylabel('Count')
		plt.xlabel('Dist. missed variants and closest SV in 1kg-ONT set')
		pdf.savefig()
		plt.close()


		plt.figure()
		print('number of matches (distance=0): ' + str(dataset.count(0)))
		print('min distance:' + str(min(data)))
		print('max distance: ' + str(max(data)))
		plt.xlim([0,5000])
		plt.hist(dataset, bins=np.arange(0, 5000 + 100, 100), alpha=0.7)
		plt.yscale('log')
		plt.ylabel('Count')
		plt.xlabel('Dist. missed variants and closest SV in 1kg-ONT set')
		pdf.savefig()
		plt.close()


		plt.figure()
		print('number of matches (distance=0): ' + str(dataset.count(0)))
		print('min distance:' + str(min(data)))
		print('max distance: ' + str(max(data)))
		plt.xlim([0,50000])
		plt.hist(dataset, bins=np.arange(0, 50000 + 100, 100), alpha=0.7)
		plt.yscale('log')
		plt.ylabel('Count')
		plt.xlabel('Dist. missed variants and closest SV in 1kg-ONT set')
		pdf.savefig()
		plt.close()




if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='bedtools closest -d -a <augmented-graph-bubbles> -b <original-graph-bubbles> | python3 plot-closest.py', description=__doc__)
	parser.add_argument('-outname', metavar='OUTNAME', default = 'histograms.pdf', help='name of the output file.')
	args = parser.parse_args()

	data = parse_tsv()
	plot_histogram(data, args.outname)	

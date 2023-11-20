import sys
import argparse
import numpy as np
from matplotlib import pyplot as plt


def parse_tsv(filename):
	segment_to_bubble = {}
	bubble_to_alignments = {}
	bubble_id = 0
	for line in open(filename, 'r'):
		fields = line.split()
		segments = [s for s in fields[11].split(',')]
		for s in segments:
			segment_to_bubble[s] = bubble_id
		bubble_id += 1
		bubble_to_alignments[bubble_id] = [(fields[0], int(fields[1]), int(fields[2])) , set([])]
	return segment_to_bubble, bubble_to_alignments

def parse_gaf(tsvs, gafs):

	for tsv, gaf in zip(tsvs, gafs):
		segment_to_bubble, bubble_to_alignments = parse_tsv(tsv)
		for line in open(filename, 'r'):
			fields = line.split()
			segments = [s for s in fields[5].split('>', '<')]
			for s in segments:
				b = segment_to_bubble[s]
				bubble_to_alignments[b][1].add(fields[0])
		print(bubble_to_alignment)

		# TODO: print number of alignments along the chromosome. Use size of bubble as bin and hight corresponds
		# to the number of alignments in that interval
		
		
	
		
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='plot-bubble-stats.py', description=__doc__)
	parser.add_argument('-files', metavar='FILE', nargs='+', required=True, help='File in TSV format.')
	parser.add_argument('-labels', metavar='LABELS', nargs='+', required=True, help='Labels. One per input file.')
	parser.add_argument('-outname', metavar='OUTNAME', default = 'histograms.pdf', help='name of the output file.')
	args = parser.parse_args()
	
	filenames = [f for f in args.files]
	labels = [l for l in args.labels]
	
	print(filenames, labels)
	if len(filenames) != len(labels):
		raise RuntimeError('Number of labels does not match number of files given.')
	
	run_plotting(filenames, labels, args.outname)

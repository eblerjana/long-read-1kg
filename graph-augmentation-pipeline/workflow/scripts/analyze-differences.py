import sys
import argparse
from collections import defaultdict

def get_reads(gaf_name):
	"""
	Extract all readnames and return a set
	that contains all of them.
	"""
	names = set([])
	for line in open(gaf_name, 'r'):
		if line.startswith("@"):
			continue
		fields = line.strip().split()
		rname = fields[0]
		names.add(rname)
	return names


def get_ranks(readname, readname_to_nodes, node_to_rank):
	result = []
	assert readname in readname_to_nodes
	nodes = readname_to_nodes[readname]
	for node in nodes:
		assert node in node_to_rank
		result.append(node_to_rank[node])
	result = set(result)
	return ','.join(list(result))


def get_coordinates(rnames1, rnames2, outname1, outname2, gfa1, gfa2, gaf1, gaf2):
	"""
	Given lists of readnames, extract their
	coordinates from input BED.
	"""
	processed = 0
	readname_to_nodes1, node_to_rank1 = parse_graph(gfa1, gaf1, rnames1)
	readname_to_nodes2, node_to_rank2 = parse_graph(gfa2, gaf2, rnames2)
	with open(outname1, 'w') as outfile1, open(outname2, 'w') as outfile2:
		for line in sys.stdin:
			if processed % 100000 == 0:
				print('Processed ' + str(processed) + ' alignments.')
			fields = line.strip().split()
			name = fields[3]
			if name in rnames1:
				outfile1.write('\t'.join([fields[0], fields[1], fields[2], fields[3], get_ranks(name, readname_to_nodes1, node_to_rank1)]) + '\n')
			if name in rnames2:
				outfile2.write('\t'.join([fields[0], fields[1], fields[2], fields[3], get_ranks(name, readname_to_nodes2, node_to_rank2)]) + '\n')
			processed += 1


def run_analysis(gaf1, gaf2, name1, name2, outname, gfa1, gfa2):

	rnames1 = get_reads(gaf1)
	print('Read ' + str(len(rnames1)) + ' reads from ' + gaf1 + '.')
	rnames2 = get_reads(gaf2)
	print('Read ' + str(len(rnames2)) + ' reads from ' + gaf2 + '.')

	outname1 = outname + '_only_' + name1 + '.bed'
	outname2 = outname + '_only_' + name2 + '.bed'

	only_gaf1 = rnames1.difference(rnames2)
	print(str(len(only_gaf1)) + ' reads are only in ' + gaf1 + '.')
	only_gaf2 = rnames2.difference(rnames1)
	print(str(len(only_gaf2)) + ' reads are only in ' + gaf2 + '.')

	get_coordinates(only_gaf1, only_gaf2, outname1, outname2, gfa1, gfa2, gaf1, gaf2)
	print('Wrote coordinates to ' + outname1 + ' and ' + outname2 + '.')


def parse_graph(gfa_file, gaf_file, readnames):
	readname_to_nodes = defaultdict(list)
	node_to_rank = {}
	for alignment in open(gaf_file, 'r'):
		fields = alignment.strip().split()
		if not fields[0] in readnames:
			continue
		node_seq = fields[5].strip().replace("<", ">")
		nodes = node_seq.split('>')
		for node in nodes:
			node_to_rank[node] = None
		readname_to_nodes[fields[0]] += [n for n in nodes if n != ""]

	for line in open(gfa_file, 'r'):
		if not line.startswith('S'):
			continue
		fields = line.strip().split()
		seq_name = ""
		seq_offset = ""
		if fields[1] not in node_to_rank:
			continue
		for field in fields[3:]:
			if field.startswith('SN:Z:'):
				seq_name = field.split(':')[-1]
			if field.startswith('SO:i:'):
				seq_offset = field.split(':')[-1]
		if seq_name != '':
			node_to_rank[fields[1]] = seq_name + ':' + seq_offset
	print(readname_to_nodes)
	print(node_to_rank)
	return readname_to_nodes, node_to_rank	
	
	



if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='bedtools bamtobed -i <BAM/CRAM> | python3 analyze-alignments.py', description=__doc__)
	parser.add_argument('gaf1', metavar='GAF1', help="First GAF file to compare.")
	parser.add_argument('gaf2', metavar='GAF2', help="Second GAF file to compare.")
	parser.add_argument('gfa1', metavar='GFA1', help="First GFA file.")
	parser.add_argument('gfa2', metavar='GFA2', help="Second GFA file.")
	parser.add_argument('-name1', metavar='NAME1', default="gaf1", help="name of first GAF.")
	parser.add_argument('-name2', metavar='NAME2', default="gaf2", help="name of second GAF.")
	parser.add_argument('-outname', metavar='OUTNAME', required=True, help="prefix of output files.")
	args = parser.parse_args()

	run_analysis(args.gaf1, args.gaf2, args.name1, args.name2, args.outname, args.gfa1, args.gfa2)

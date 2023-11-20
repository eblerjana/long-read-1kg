import sys

#ref = sys.argv[1]

#for line in sys.stdin:
#	sample = line.strip()
#	line = [sample, "/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/" + sample + "/pav_" + ref + "/pav_svtigs_merged.vcf"]
#	print("\t".join(line))


# /gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/HG03521/pav_hg38/pav_svtigs_merged.vcf

for line in sys.stdin:
	fields = line.split('/')
	sample = fields[-3]
	print('\t'.join([sample, line.strip()]))

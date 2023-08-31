configfile: "config/config.yaml"
callsets = [c for c in config["callset_vcfs"].keys()]
chromosomes = [i for i in config["minigraph_gfa"].keys()]


checkpoint create_paths:
	"""
	Given a multisample callset VCF, split variants in non-overlapping subsets.
	Produces a new VCF per subset with haploid haplotype path.
	"""
	input:
		lambda wildcards: config['callset_vcfs'][wildcards.callset]
	output:
		directory("results/paths/{callset}/")
	resources:
		mem_total_mb=10000,
		runtime_hrs=1
	wildcard_constraints:
		callset = "|".join(callsets)
	log:
		"results/paths/{callset}/{callset}.log"
	benchmark:
		"results/paths/{callset}/{callset}-benchmark.txt"
	params:
		outprefix = "results/paths/{callset}/{callset}"
	shell:
		"""
		python3 workflow/scripts/create_paths.py -vcf {input} -single {params.outprefix} &> {log}
		"""



rule compress_vcf:
	input:
		"results/paths/{callset}/{callset}_path{path_id}.vcf"
	output:
		"results/paths/{callset}/{callset}_path{path_id}.vcf.gz"
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		bgzip -c {input} > {output}
		tabix -p vcf {output}
		"""




rule extract_region_fasta:
	input:
		config['reference']
	output:
		fasta=temp("results/reference/reference_{chrom}.fa"),
		fai=temp("results/reference/reference_{chrom}.fa.fai")
	conda:
		"../envs/minigraph.yml"
	log:
		"results/reference/reference_{chrom}.log"
	benchmark:
		"results/reference/reference_{chrom}_benchmark.txt"
	shell:
		"""
		samtools faidx {input} {wildcards.chrom} > {output.fasta}
		samtools faidx {output.fasta}
		"""



rule compute_consensus_region:
	"""
	Given the haploid VCFs produced for each path, construct
	a consensus sequence by inserting them into the reference
	genome for specified region
	"""
	input:
		reference="results/reference/reference_{chrom}.fa",
		vcf = "results/paths/{callset}/{callset}_path{path_id}.vcf.gz"
	output:
		fasta = temp("results/paths/{callset}/{callset}_path{path_id}_{chrom}.fa"),
		tmp = temp("results/paths/{callset}/tmp/{callset}_path{path_id}_{chrom}.vcf.gz"),
		tbi = temp("results/paths/{callset}/tmp/{callset}_path{path_id}_{chrom}.vcf.gz.tbi")
	log:
		"results/paths/{callset}/{callset}_path{path_id}_{chrom}_consensus.log"
	wildcard_constraints:
		callset = "|".join(callsets),
		chrom = "|".join([c for c in chromosomes if c != "all"])
	benchmark:
		"results/paths/{callset}/{callset}_path{path_id}_{chrom}_consensus_benchmark.txt"
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		bcftools view -r {wildcards.chrom} {input.vcf} | bgzip > {output.tmp}
		tabix -p vcf {output.tmp}
		bcftools consensus -f {input.reference} {output.tmp} 2> {log} | python3 workflow/scripts/rename-fasta.py {wildcards.callset}_{wildcards.path_id} > {output.fasta}
		"""



rule compute_consensus_full:
	"""
	Given the haploid VCFs produced for each path, construct
	a consensus sequence by inserting them into the reference
	genome.
	"""
	input:
		reference=config['reference'],
		vcf = "results/paths/{callset}/{callset}_path{path_id}.vcf.gz"
	output:
		fasta = temp("results/paths/{callset}/{callset}_path{path_id}_{chrom}.fa")
	log:
		"results/paths/{callset}_path{path_id}_{chrom}_consensus.log"
	wildcard_constraints:
		callset = "|".join(callsets),
		chrom = "all"
	benchmark:
		"results/paths/{callset}/{callset}_path{path_id}_{chrom}_consensus_benchmark.txt"
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		bcftools consensus -f {input.reference} {input.vcf} 2> {log} | python3 workflow/scripts/rename-fasta.py {wildcards.callset}_{wildcards.path_id} > {output.fasta}
		"""



ruleorder: create_paths > extend_minigraph

def aggregate_input(wildcards):
	outfiles = []
	for callset in callsets:
		checkpoint_output = checkpoints.create_paths.get(callset=callset).output[0]
		outfiles += [i for i in  expand("results/paths/" + callset +  "/" + callset + "_path{path_id}_{chrom}.fa", chrom=wildcards.chrom, path_id = glob_wildcards(os.path.join(checkpoint_output, callset + "_path{path_id}.vcf")).path_id)]
	return outfiles

rule extend_minigraph:
	"""
	Add each new consensus sequence to the minigraph.
	"""
	input:
		paths = aggregate_input,
		minigraph = lambda wildcards: config["minigraph_gfa"][wildcards.chrom]
	output:
		"results/minigraph/minigraph-extended_{chrom}.gfa"
	log:
		"results/minigraph/minigraph_{chrom}.log"
	resources:
		mem_total_mb=100000,
		runtime_hrs=10
	benchmark:
		"results/minigraph/minigraph_{chrom}_benchmark.txt"
	threads: 24
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		minigraph -cxggs -t{threads} {input.minigraph} {input.paths} 2> {log} 1> {output} 
		"""


rule gfa_bubble_stats:
	"""
	Compute number of bubbles in the GFA files before/after inserting variants
	"""
	input:
		lambda wildcards: config["minigraph_gfa"][wildcards.chrom] if wildcards.version == "original" else "results/minigraph/minigraph-extended_{chrom}.gfa"
	output:
		"results/statistics/{version}-graph-{chrom}_bubbles.tsv"
	conda:
		"../envs/minigraph.yml"
	wildcard_constraints:
		version = "original|extended"
	shell:
		"gfatools bubble {input} > {output}"

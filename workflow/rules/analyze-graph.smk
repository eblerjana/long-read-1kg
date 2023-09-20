configfile: "config/config.yaml"
chromosomes = [i for i in config["minigraph_gfa"].keys()] 
gaftools = config['gaftools']


#############################################
#   Count bubble structures in the graph
#############################################


rule gfa_bubble_stats:
	"""
	Compute number of bubbles in the GFA files before/after inserting variants
	"""
	input:
		lambda wildcards: config["minigraph_gfa"][wildcards.chrom] if wildcards.version == "original" else "results/minigraph/minigraph-extended_{chrom}.gfa"
	output:
		"results/statistics/bubbles/{version}_{chrom}_bubbles.tsv"
	conda:
		"../envs/minigraph.yml"
	wildcard_constraints:
		version = "original|extended"
	resources:
		mem_total_mb=50000
	shell:
		"gfatools bubble {input} > {output}"


rule gaftools_order_gfa:
	"""
	Compute bubble statistics with gaftools order_gfa
	"""
	input:
		lambda wildcards: config["minigraph_gfa"][wildcards.chrom] if wildcards.version == "original" else "results/minigraph/minigraph-extended_{chrom}.gfa"
	output:
		directory("results/statistics/ordering/{version}_{chrom}/")
	log:
		"results/statistics/{version}/ordering/{version}_{chrom}_ordering.log"
	wildcard_constraints:
		version = "original|extended"
	shell:
		"""
		{gaftools} order_gfa --outdir {output} {input} &> {log}
		"""



#############################################
#   Check read mappability
#############################################


rule concat_gfas:
	input:
		lambda wildcards: [config["minigraph_gfa"][c] for c in chromosomes] if wildcards.version == "original" else ["results/minigraph/minigraph-extended_" + c + ".gfa" for c in chromosomes]
	output:
		temp("results/minigraph/full/minigraph-{version}-full.gfa")
	shell:
		"""
		cat {input} > {output}
		"""


rule minigraph_align:
	"""
	Align reads with minigraph to the graphs
	"""
	input:
		graph = "results/minigraph/full/minigraph-{version}-full.gfa",
		reads = lambda wildcards: config["reads"][wildcards.sample]
	output:
		"results/statistics/mapping/{version}_all_{sample}_minigraph.gaf"
	log:
		"results/statistics/mapping/{version}_all_{sample}_minigraph.log"
	benchmark:
		"results/statistics/mapping/{version}_all_{sample}_minigraph_benchmark.txt"
	wildcard_constraints:
		version = "original|extended"
	resources:
		mem_total_mb=200000,
		runtime_hrs=5
	threads: 24
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		minigraph -cx lr -t{threads} {input.graph} {input.reads} 2> {log} 1> {output}
		"""


rule graphaligner_align:
	"""
	Align reads with Graphaligner to the graphs
	"""
	input:
		graph = "results/minigraph/full/minigraph-{version}-full.gfa",
		reads = lambda wildcards: config["reads"][wildcards.sample]
	output:
		"results/statistics/mapping/{version}_all_{sample}_graphaligner.gaf"
	log:
		"results/statistics/mapping/{version}_all_{sample}_graphaligner.log"
	benchmark:
		"results/statistics/mapping/{version}_all_{sample}_graphaligner_benchmark.txt"
	wildcard_constraints:
		version = "original|extended"
	resources:
		mem_total_mb=200000,
		runtime_hrs=10
	threads: 24
	conda:
		"../envs/graphaligner.yml"
	shell:
		"""
		GraphAligner -g {input.graph} -f {input.reads} -a {output} -t {threads} -x vg &> {log}
		"""




rule alignment_stats:
	"""
	Compute alignment statistics.
	"""
	input:
		"results/statistics/mapping/{version}_all_{sample}_{method}.gaf"
	output:
		"results/statistics/mapping/{version}_all_{sample}_{method}.stats"
	wildcard_constraints:
		version = "original|extended",
		method = "minigraph|graphaligner"
	resources:
		mem_total_mb=10000
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		{gaftools} stat {input} &> {output}
		"""

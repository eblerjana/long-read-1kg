configfile: "config/config.yaml"
gaftools = config['gaftools']


#############################################
#   Count bubble structures in the graph
#############################################


rule gfa_bubble_stats:
	"""
	Compute number of bubbles in the GFA files before/after inserting variants
	"""
	input:
		lambda wildcards: config["minigraph_gfa"] if wildcards.version == "original" else "results/minigraph/minigraph-extended_all.gfa"
	output:
		"results/statistics/bubbles/{version}_all_bubbles.tsv"
	conda:
		"../envs/minigraph.yml"
	wildcard_constraints:
		version = "original|extended"
	resources:
		mem_total_mb=50000
	shell:
		"gfatools bubble {input} > {output}"


rule plot_bubble_sizes:
	"""
	Plot histogram of bubble sizes
	"""
	input:
		original = "results/statistics/bubbles/original_all_bubbles.tsv",
		extended = "results/statistics/bubbles/extended_all_bubbles.tsv"
	output:
		"results/statistics/plots/all_bubbles.pdf"
	shell:
		"python3 workflow/scripts/plot-histogram.py -files {input.original} {input.extended} -labels original_graph_all augmented_graph_all -outname {output}"



rule gaftools_order_gfa:
	"""
	Compute bubble statistics with gaftools order_gfa
	"""
	input:
		lambda wildcards: config["minigraph_gfa"] if wildcards.version == "original" else "results/minigraph/minigraph-extended_all.gfa"
	output:
		directory("results/statistics/ordering/{version}_all/")
	log:
		"results/statistics/{version}/ordering/{version}_all_ordering.log"
	wildcard_constraints:
		version = "original|extended"
	shell:
		"""
		{gaftools} order_gfa --outdir {output} {input} &> {log}
		"""



#############################################
#   Check read mappability
#############################################


rule minigraph_align:
	"""
	Align reads with minigraph to the graphs
	"""
	input:
		graph = lambda wildcards: config["minigraph_gfa"] if wildcards.version == "original" else "results/minigraph/minigraph-extended_all.gfa",
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
		minigraph -cx lr -t{threads} {input.graph} {input.reads} --vc 2> {log} 1> {output}
		"""


rule graphaligner_align:
	"""
	Align reads with Graphaligner to the graphs
	"""
	input:
		graph = lambda wildcards: config["minigraph_gfa"] if wildcards.version == "original" else "results/minigraph/minigraph-extended_all.gfa",
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



rule alignment_read_stats:
	"""
	Compute read-alignment statistics.
	"""
	input:
		"results/statistics/mapping/{version}_all_{sample}_{method}.gaf"
	output:
		"results/statistics/mapping/{version}_all_{sample}_{method}_reads.stats"
	wildcard_constraints:
		version = "original|extended",
		method = "minigraph|graphaligner"
	resources:
		mem_total_mb=10000
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		cat {input} | python3 workflow/scripts/analyze-alignments.py &> {output}
		"""

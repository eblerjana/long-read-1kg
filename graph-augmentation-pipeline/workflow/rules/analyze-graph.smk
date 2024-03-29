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


rule plot_distance_histogram:
	"""
	Plot histogram of bubble distances
	"""
	input:
		original = "results/statistics/bubbles/original_all_bubbles.tsv",
		extended = "results/statistics/bubbles/extended_all_bubbles.tsv"
	output:
		tsv="results/statistics/plots/bubble_distances.tsv",
		pdf="results/statistics/plots/bubble_distances.pdf"
	log:
		"results/statistics/plots/bubble_distances.log"
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		bedtools closest -d -t first -a {input.extended} -b {input.original} > {output.tsv}
		cat {output.tsv} | python3 workflow/scripts/plot-distances.py -outname {output.pdf} &> {log}
		"""


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


rule extract_raw_reads:
	"""
	From the given bam/cram file, extract raw reads.
	"""
	input:
		alignments = lambda wildcards: config["reads"][wildcards.sample],
		reference = config["reference"]
	output:
		temp("results/statistics/mapping/raw-reads/{sample}_raw.fasta")
	conda:
		"../envs/minigraph.yml"
	resources:
		mem_total_mb = 50000
	shell:
		"""
		samtools fasta --reference {input.reference} {input.alignments} > {output}
		"""


rule minigraph_align:
	"""
	Align reads with minigraph to the graphs
	"""
	input:
		graph = lambda wildcards: config["minigraph_gfa"] if wildcards.version == "original" else "results/minigraph/minigraph-extended_all.gfa",
		reads = "results/statistics/mapping/raw-reads/{sample}_raw.fasta"
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
		reads = "results/statistics/mapping/raw-reads/{sample}_raw.fasta"
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


rule graph_stats:
	"""
	Compute graph statistics.
	"""
	input:
		graph = lambda wildcards: config["minigraph_gfa"] if wildcards.version == "original" else "results/minigraph/minigraph-extended_all.gfa"
	output:
		"results/statistics/graph/{version}_graph-stats.txt"
	wildcard_constraints:
		version = "original|extended"
	shell:
		"""
		cat {input} | python3 workflow/scripts/analyze-graph.py &> {output}
		"""


rule analyze_alignment_differences:
	"""
	Analyze reads that only align in original/extended
	graph.
	"""
	input:
		original_gaf = "results/statistics/mapping/original_all_{sample}_{method}.gaf",
		extended_gaf = "results/statistics/mapping/extended_all_{sample}_{method}.gaf",
		original_gfa = config["minigraph_gfa"],
		extended_gfa = "results/minigraph/minigraph-extended_all.gfa",
		reads = lambda wildcards: config["reads"][wildcards.sample],
		reference = config["reference"]
	output:
		original = "results/statistics/mapping/alignments_{sample}_{method}_only_original.bed",
		extended = "results/statistics/mapping/alignments_{sample}_{method}_only_extended.bed"
	log:
		"results/statistics/mapping/alignments_{sample}_{method}.log"
	wildcard_constraints:
		method = "minigraph|graphaligner"
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 2
	conda:
		"../envs/minigraph.yml"
	params:
		outname = "results/statistics/mapping/alignments_{sample}_{method}"
	shell:
		"""
		samtools view -h --reference {input.reference} {input.reads} | bedtools bamtobed -i -  | python3 workflow/scripts/analyze-differences.py {input.original_gaf} {input.extended_gaf} {input.original_gfa} {input.extended_gfa} -name1 original -name2 extended -outname {params.outname} &> {log}
		"""

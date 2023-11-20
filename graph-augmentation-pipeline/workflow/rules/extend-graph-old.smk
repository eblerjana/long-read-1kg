configfile: "config/config.yaml"
callsets = [c for c in config["callset_vcfs"].keys()]


rule create_paths:
	input:
		lambda wildcards: config['callset_vcfs'][wildcards.callset]
	output:
		dynamic("results/paths/{callset}/{callset}_path{path_id}.vcf.gz"),
	resources:
		mem_total_mb=10000,
		runtime_hrs=1
	wildcard_constraints:
		callset = "|".join(callsets)
	log:
		"results/paths/{callset}.log"
	params:
		outprefix = "results/paths/{callset}"
	shell:
		"""
		python3 workflow/scripts/create_paths.py -vcf {input} -single {params.outprefix} &> {log}
		"""


rule compute_consensus:
	input:
		reference=config['reference'],
		vcf = "results/paths/{callset}/{callset}_path{path_id}.vcf.gz"
	output:
		temp("results/paths/{callset}/{callset}_path{path_id}.fa")
	log:
		"results/paths/{callset}_path{path_id}_consensus.log"
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		bcftools consensus -f {input.reference} {input.vcf} 2> {log} 1> {output}
		"""

rule extend_minigraph:
	input:
		paths = dynamic(expand("results/paths/{callset}/{callset}_path{{path_id}}.fa", callset = callsets)),
		minigraph = config["minigraph_gfa"]
	output:
		"results/minigraph/minigraph-extended.gfa"
	log:
		"results/minigraph/minigraph.log"
	resources:
		mem_total_mb=50000,
		runtime_hrs=5
	threads: 24
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		minigraph -cxggs -t{threads} {input.minigraph} {input.paths} 2> {log} 1> {output} 
		"""
	





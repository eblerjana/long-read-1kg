configfile: "config/config.yaml"
callsets = [c for c in config["callset_vcfs"].keys()]


checkpoint create_paths:
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
		bcftools consensus -f {input.reference} {input.vcf} 2> {log} | python3 workflow/scripts/rename-fasta.py {wildcards.callset}_{wildcards.path_id} > {output}
		"""

ruleorder: create_paths > extend_minigraph


def aggregate_input(wildcards):
	outfiles = []
	for callset in callsets:
		checkpoint_output = checkpoints.create_paths.get(callset=callset).output[0]
		outfiles += [i for i in  expand("results/paths/" + callset +  "/" + callset + "_path{path_id}.fa", path_id = glob_wildcards(os.path.join(checkpoint_output, callset + "_path{path_id}.vcf")).path_id)]
	print("outfiles", outfiles)
	return outfiles


rule extend_minigraph:
	input:
		paths = aggregate_input,
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
	





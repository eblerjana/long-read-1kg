configfile: "config/config.yaml"
callsets = [c for c in config["callset_vcfs"].keys()]

margin = 500 # (= 1000)


rule remove_centromer_variants:
	"""
	Remove variants in centromer regions from the VCF.
	"""
	input:
		vcf = lambda wildcards: config['callset_vcfs'][wildcards.callset],
		bed = config["mask"]
	output:
		"results/masked_vcfs/{callset}_masked.vcf.gz"
	log:
		"results/masked_vcfs/{callset}_masked.log"
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		bedtools subtract -a {input.vcf} -b {input.bed} -A 2> {log} | bgzip > {output}
		tabix -p vcf {output}
		"""

rule mask_reference:
	"""
	Set centromer regions to N.
	"""
	input:
		fasta = config["reference"],
		bed = config["mask"]
	output:
		"results/masked_reference/reference_masked.fa"
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output}
		samtools faidx {output}
		"""


checkpoint create_paths:
	"""
	Given a multisample callset VCF, split variants in non-overlapping subsets.
	Produces a new VCF per subset with haploid haplotype path.
	"""
	input:
		"results/masked_vcfs/{callset}_masked.vcf.gz"
	output:
		directory("results/paths/{callset}_vcfs/")
	resources:
		mem_total_mb=10000,
		runtime_hrs=1
	wildcard_constraints:
		callset = "|".join(callsets)
	log:
		"results/paths/{callset}_vcfs/{callset}.log"
	benchmark:
		"results/paths/{callset}_vcfs/{callset}-benchmark.txt"
	params:
		outprefix = "results/paths/{callset}_vcfs/{callset}"
	shell:
		"""
		python3 workflow/scripts/create_paths.py -vcf {input} -single {params.outprefix} -margin {margin} &> {log}
		"""



rule compress_vcf:
	input:
		"results/paths/{callset}_vcfs/{callset}_path{path_id}.vcf"
	output:
		"results/paths/{callset}_gz/{callset}_path{path_id}.vcf.gz"
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		bgzip -c {input} > {output}
		tabix -p vcf {output}
		"""


rule compute_consensus_full:
	"""
	Given the haploid VCFs produced for each path, construct
	a consensus sequence by inserting them into the reference
	genome. Rename FASTA records to ensure that each FASTA file
	has unique names (required by minigraph). Finally, remove all
	Ns from sequences and split FASTA record at each N.
	"""
	input:
		reference="results/masked_reference/reference_masked.fa",
		vcf = "results/paths/{callset}_gz/{callset}_path{path_id}.vcf.gz",
	output:
		fasta_tmp1 = temp("results/paths/{callset}_fasta/tmp1_{callset}_path{path_id}_all.fa"),
		fasta = "results/paths/{callset}_fasta/{callset}_path{path_id}_all.fa"
	log:
		cons = "results/paths/{callset}_fasta/{callset}_path{path_id}_all_consensus.log",
		split = "results/paths/{callset}_fasta/{callset}_path{path_id}_all_split.log"
	wildcard_constraints:
		callset = "|".join(callsets)
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 4
	benchmark:
		"results/paths/{callset}_fasta/{callset}_path{path_id}_all_consensus_benchmark.txt"
	conda:
		"../envs/minigraph.yml"
	params:
		name = "{callset}_{path_id}_all"
	shell:
		"""
			bcftools consensus -f {input.reference} {input.vcf} 2> {log.cons} | python3 workflow/scripts/rename-fasta.py {wildcards.callset}_{wildcards.path_id} > {output.fasta_tmp1}
			cat {output.fasta_tmp1} | python3 workflow/scripts/split_fasta.py -o {output.fasta} > {log.split}
		"""



ruleorder: create_paths > extend_minigraph

def aggregate_input(wildcards):
	outfiles = []
	for callset in callsets:
		checkpoint_output = checkpoints.create_paths.get(callset=callset).output[0]
		outfiles += [i for i in  expand("results/paths/" + callset +  "_fasta/" + callset + "_path{path_id}_all.fa", path_id = glob_wildcards(os.path.join(checkpoint_output, callset + "_path{path_id}.vcf")).path_id)]
	return sorted(outfiles)

rule extend_minigraph:
	"""
	Add each new consensus sequence to the minigraph.
	"""
	input:
		paths = aggregate_input,
		minigraph = config["minigraph_gfa"]
	output:
		"results/minigraph/minigraph-extended_all.gfa"
	log:
		"results/minigraph/minigraph_all.log"
	resources:
		mem_total_mb=200000,
		runtime_hrs=30
	benchmark:
		"results/minigraph/minigraph_all_benchmark.txt"
	threads: 32
	conda:
		"../envs/minigraph.yml"
	shell:
		"""
		minigraph -cxggs -t{threads} {input.minigraph} {input.paths} 2> {log} 1> {output} 
		"""

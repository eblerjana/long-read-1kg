configfile: "config/config.yaml"

callsets = [c for c in config["callset_vcfs"].keys()]
max_af = config['max_af']
min_af = config['min_af']
normalize = config["reference"] != ""


rule normalize_vcf:
	input:
		vcf = lambda wildcards: config["callset_vcfs"][wildcards.callset],
		reference = config["reference"]
	output:
		"results/vcfs/{callset}-normalized.vcf.gz"
	conda:
		"../envs/comparison.yml"
	wildcard_constraints:
		callset = "|".join(callsets)
	shell:
		"""
		bcftools norm -c x -f {input.reference} {input.vcf} -Oz -o {output}
		"""


rule vcf_stats:
	input:
		lambda wildcards: "results/vcfs/{callset}-normalized.vcf.gz" if normalize else config["callset_vcfs"][wildcards.callset]
	output:
		"results/vcf-stats/{callset}.stats"
	conda:
		"../envs/comparison.yml"
	wildcard_constraints:
		callset = "|".join(callsets)
	shell:
		"zcat {input} | python3 workflow/scripts/set-pass.py | bcftools view -f PASS --min-af {min_af} --max-af {max_af} | python3 workflow/scripts/vcf_stats.py > {output}"



rule add_tags:
	input:
		lambda wildcards: "results/vcfs/{callset}-normalized.vcf.gz" if normalize else config["callset_vcfs"][wildcards.callset]
	output:
		"results/vcfs/{callset}-tagged.vcf"
	wildcard_constraints:
		callset = "|".join(callsets)
	log:
		"results/vcfs/{callset}-tagged.log"
	conda:
		"../envs/comparison.yml"
	params:
		ignore_ids = "--ignore-ids" if normalize else ""
	shell:
		"zcat {input} | python3 workflow/scripts/set-pass.py | bcftools view -f PASS --min-af {min_af} --max-af {max_af} | python3 workflow/scripts/add-svtags.py {params.ignore_ids} 2> {log} 1> {output}"


rule extract_sample:
	input:
		"results/vcfs/{callset}-tagged.vcf"
	output:
		"results/vcfs/{callset}-tagged-{sample}.vcf"
	conda:
		"../envs/comparison.yml"
	shell:
		"bcftools view --samples {wildcards.sample} {input} | bcftools view --min-ac 1 > {output}"


def intersect_vcfs_files(wildcards):
	files = []
	for c in COMBINATIONS[wildcards.combination]:
		if wildcards.sample == "all":
			files.append("results/vcfs/{callset}-tagged.vcf".format(callset = c))
		else:
			files.append("results/vcfs/{callset}-tagged-{sample}.vcf".format(callset = c, sample = wildcards.sample))
	return files


rule intersect_vcfs:
	input:
		intersect_vcfs_files
	output:
		tsv="results/{combination}/intersection/intersection_{combination}_{sample}.tsv",
		vcf="results/{combination}/intersection/intersection_{combination}_{sample}.vcf",
		pdf="results/{combination}/intersection/intersection_{combination}_{sample}.pdf",
		plot="results/{combination}/intersection/intersection_upset_{combination}_{sample}.pdf"
	conda:
		"../envs/upsetplot.yml"
	log:
		intersect="results/{combination}/intersection/intersection_{combination}_{sample}.log",
		plot="results/{combination}/intersection/plotting_{combination}_{sample}.log"
	params:
		names =  lambda wildcards: [c for c in COMBINATIONS[wildcards.combination]],
		columns = lambda wildcards: ["in_" + c for c in COMBINATIONS[wildcards.combination]]
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 2
	shell:
		"""
		python3 workflow/scripts/intersect_callsets.py intersect -c {input} -n {params.names} -t {output.tsv} -v {output.vcf} -p {output.pdf} --id-from-vcf &> {log.intersect}
		python3 workflow/scripts/plot-upset.py -t {output.tsv} -o {output.plot} -n {params.columns} &> {log.plot}
		"""

rule combined_plot:
	input:
		expand("results/{{combination}}/intersection/intersection_{{combination}}_{sample}.tsv", sample = SAMPLES)
	output:
		pdf = "results/{combination}/intersection/intersection_bar_{combination}.pdf",
		tsv = "results/{combination}/intersection/intersection_summary_{combination}.tsv"
	conda:
		"../envs/upsetplot.yml"
	shell:
		"""
		ls {input} | python3 workflow/scripts/plot-intersections.py -o {output.pdf} &> {output.tsv}
		"""	







rule get_all_variants:
	input:
		lambda wildcards: INTERSECTIONS[wildcards.sample]['table']
	output:
		"results/{sample}/{sample}_all-variants.bed"
	params:
		present = lambda wildcards: INTERSECTIONS[wildcards.sample]['absent']
	shell:
		"""
		cat {input} | python3 workflow/scripts/get_subset.py --present {params.present} > {output} 
		"""


rule get_missed_variants:
	input:
		lambda wildcards: INTERSECTIONS[wildcards.sample]['table']
	output:
		"results/{sample}/{sample}_missed-variants.bed"
	params:
		present = lambda wildcards: INTERSECTIONS[wildcards.sample]['present'],
		absent = lambda wildcards: INTERSECTIONS[wildcards.sample]['absent']
	shell:
		"""
		cat {input} | python3 workflow/scripts/get_subset.py --present {params.present} --absent {params.absent} > {output}
		"""


rule prepare_panel_regions:
	input:
		lambda wildcards: PANEL_VCF if wildcards.region == 'panel' else PANEL_BI_VCF
	output:
		"results/{region}-regions.bed"
	wildcard_constraints:
		region = "panel|sites"
	shell:
		"zcat {input} | python3 workflow/scripts/vcf-to-bed.py | bedtools merge > {output}"


rule prepare_graph_regions:
	input:
		GRAPH_BUBBLES
	output:
		"results/graph-regions.bed"
	shell:
		"""
		cut -f1,2,3 {input} > {output}
		"""

rule prepare_genotypes:
	input:
		lambda wildcards: UNFILTERED if wildcards.region == "unfiltered-gt" else FILTERED
	output:
		"results/{region}-regions.bed"
	wildcard_constraints:
		region = "filtered-gt|unfiltered-gt"
	shell:
		"zcat {input} | python3 workflow/scripts/vcf-to-bed.py | bedtools merge > {output}" 


rule compare_to_regions:
	input:
		variants = "results/{sample}/{sample}_missed-variants.bed",
		regions = "results/{region}-regions.bed"
	output:
		inside="results/{sample}/{sample}_missed-inside-{region}.bed",
		outside="results/{sample}/{sample}_missed-outside-{region}.bed"
	shell:
		"""
		bedtools intersect -wa -u -a {input.variants} -b {input.regions} > {output.inside}
		bedtools subtract -A -a {input.variants} -b {input.regions} > {output.outside}
		"""


rule compute_closest:
	input:
		variants = "results/{sample}/{sample}_missed-variants.bed",
		all = "results/{sample}/{sample}_all-variants.bed"
	output:
		"results/{sample}/{sample}_missed-variants-closest.bed"
	shell:
		"""
		bedtools closest -a {input.variants} -b {input.all} -d > {output}
		"""

rule plot_closest_distance:
	input:
		"results/{sample}/{sample}_missed-variants-closest.bed"
	output:
		"results/{sample}/{sample}_missed-variants-closest.pdf"
	log:
		"results/{sample}/{sample}_missed-variants-closest.pdf.log"
	shell:
		"""
		cat {input} | python3 workflow/scripts/plot-closest.py -outname {output} &> {log}
		"""


########################## analyze overlap with genomic regions ################################


rule get_baseline_variants:
	input:
		lambda wildcards: INTERSECTIONS[wildcards.sample]['table']
	output:
		"{results}/{sample}/{sample}_baseline.bed"
	params:
		baseline = lambda wildcards: INTERSECTIONS[wildcards.sample]['baseline']
	shell:
		"""
		cat {input} | python3 workflow/scripts/get_subset.py --present {params.baseline} > {output}
		"""

rule overlap_regions:
	input:
		bed = lambda wildcards: "results/{sample}/{sample}_baseline.bed" if wildcards.set == "baseline" else "results/{sample}/{sample}_missed-variants.bed",
		region = lambda wildcards: BED[wildcards.region]
	output:
		inside = "results/{sample}/{region}/{sample}_{set}_{region}_inside.bed",
		outside = "results/{sample}/{region}/{sample}_{set}_{region}_outside.bed"
	shell:
		"""
		bedtools intersect -wa -u -a {input.bed} -b {input.region} > {output.inside}
		bedtools subtract -A -a {input.bed} -b {input.region} > {output.outside}
		"""
		

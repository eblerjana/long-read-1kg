# Graph augmentation pipeline

This pipeline extends an existing minigraph pangenome graph by adding variant calls produced by different callsets (such as Delly, Sniffles2 or CuteSV).

## Inputs

* **multi-sample VCF** files with sequence resolved variant calls across samples
* **minigraph GFA** files containing existing graph to be extended
* **reference FASTA** containing the reference sequence underlying the callsets

## Outputs

* **minigraph GFA** files with variants added

## What the pipeline does

![Pipeline overview](pipeline.png)

* **Step 1:** remove variation in centromer regions and mask centromeres in reference genome
* **Step 2:** group non-overlapping variants and prduce "pseudo" haplotypes. For each such pseudo haplotype, a VCF file is produced containing the variants it covers. Next, a consensus sequence is constructed for each such pseudo haplotype by inserting the variants into the reference genome.
* **Step 3:** Add these new consensus sequences to the GFAs using minigraph

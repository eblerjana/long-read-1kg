# How to re-produce SAGA augmented graph

The SAGA results: https://www.biorxiv.org/content/10.1101/2024.04.18.590093v1 can be re-produced with the following steps:

## Steps to re-produce

## 1. Download data

### Callsets

* sniffles:
 ``` bat
TODO
```
* delly:
 ``` bat
TODO
```
* SVarp: run pipeline: https://github.com/eblerjana/long-read-1kg/tree/main/prepare-SVarp-callset and provide the output VCF

### Minigraph
``` bat
wget -p data/ https://zenodo.org/records/6983934/files/chm13-90c.r518.gfa.gz?download=1
gunzip data/chm13-90c.r518.gfa.gz
```
### Reference

``` bat
wget -P data/ http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/reference/1KG_ONT_VIENNA_hg38.fa.gz
gunzip data/1KG_ONT_VIENNA_hg38.fa.gz
```
### Reads
``` bat
wget -P data/ http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG00513.hg38.cram
```

### Install gaftools

See https://github.com/marschall-lab/gaftools

## 2. Prepare config file

``` bat

# callsets with variants to be added to the graph
callset_vcfs:
 sniffles: "data/delins.sniffles.hg38.liftedT2T.13Nov2023.nygc.vcf.gz"
 delly: "data/delins.delly.hg38.liftedT2T.13Nov2023.nygc.vcf.gz"
 SVarp: "data/SVarp_biallelic_truvari_filtered.vcf.gz"

# minigraph GFA of full genome.
minigraph_gfa: "data/chm13-90c.r518.gfa"

# reference sequence underlying the callsets
reference: "data/1KG_ONT_VIENNA_hg38.fa"

# BED to mask reference
mask: "data/hg38-ucsc-centromers.bed"

# path to gaftools executable
gaftools: "gaftools"


#######################################################
# optional arguments (needed for evaluation only)
####################################################### 

# reads of the samples to be used for validation. Files must be in SAM/BAM/CRAM format.
# # reads must be aligned to the same reference as listed above.
reads:
 HG00513: "data/HG00513.hg38.cram"

```
## 3. Run pipeline

``` bat
snakemake --use-conda -j <number of cores>
```


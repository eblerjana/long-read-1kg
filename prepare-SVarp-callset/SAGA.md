# Steps to re-produce SAGA results

The SVarp merged calls used for graph augmentation in https://www.biorxiv.org/content/10.1101/2024.04.18.590093v1 can be re-produced using the following steps.

## 1. Download data

### SVarp individual calls (hg38)
Merged SVarp calls generated relative to GRCh38 were produced with the following pipeline: https://github.com/marschall-lab/project-ont-1kg/tree/main/snakemake_pipelines/pre-augmentation-giggles-svarp.
Generate a file (``paths.tsv``) with the following format specifying paths to the VCFs `` result/svarp-giggles/svarp/<sample>/pav_hg38/pav_svtigs_merged.vcf`` produced by the pipeline:

``` bat
<sample name>  /path/to/pav_svtigs_merged.vcf

```

### Reference

``` bat
wget -P data/ http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/reference/1KG_ONT_VIENNA_hg38.fa.gz
gunzip data/1KG_ONT_VIENNA_hg38.fa.gz
```

## 2. Prepare config

``` bat
# tsv file of format: <sample-name> <path/to/callset.vcf>
samples_list: "paths.tsv"

# reference FASTA
reference: "data/1KG_ONT_VIENNA_hg38.fa"

# how the output folder/file shall be named
results: "results"
```

## 3. Run pipeline

``` bat
snakemake --use-conda -j <number of cores>
```

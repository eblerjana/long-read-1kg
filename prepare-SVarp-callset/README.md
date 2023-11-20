# Pipeline to create multi-sample VCF from single-sample SVarp calls

Given single-sample VCFs produced from SVarp calls, this pipeline merges them into a multi-sample VCF. 


## Inputs

* a tsv file specifying the paths to the single-sample VCFs following the format:

``` bat

<sample1>    </path/to/sample1-calls.vcf.gz>
<sample2>    </path/to/sample2-calls.vcf.gz>
...

```

* reference genome underlying the calls

The input files need to be specified in the config file: `` config.yaml``. It has the following fields:

``` bat
# tsv file of format: <sample-name> <path/to/callset.vcf>
samples_list: "/path/to/callset.tsv"

# reference FASTA
reference: "/path/to/reference.fa"

# how the output folder/file shall be named
results: "results"
```

Check out the file `` config.yaml `` for an example.


## Outputs

A merged and filtered multi-sample VCF file: `` <outname>/merged-vcfs/<outname>_biallelic_truvari_filtered.vcf.gz ``


## What the pipeline does

The main steps include:

* merge files using `` bcftools merge ``
* collapse similar SV alleles using `` truvari collapse `` (https://github.com/ACEnglish/truvari/wiki/collapse)


## How to run:

Prepare the config file `` config/config.yaml `` as explained above. Then, run the pipeline using the following command:

``` bat
snakemake --use-conda -j <number of cores>
```

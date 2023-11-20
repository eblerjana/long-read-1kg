# Pipeline to create multi-sample VCF from single-sample SVarp calls

Given single-sample VCFs produced from SVarp calls, this pipeline merges them into a multi-sample VCF. Steps include:

* merge files using `` bcftools merge ``
* collapse similar SV alleles using `` truvari collapse `` (https://github.com/ACEnglish/truvari/wiki/collapse)


## Inputs

* a tsv file specifying the paths to the single-sample VCFs following the format:

``` bat

<sample1>    </path/to/sample1-calls.vcf.gz>
<sample2>    </path/to/sample2-calls.vcf.gz>
...

```

* reference genome underlying the calls


## Outputs

A merged and filtered multi-sample VCF file: `` <outname>/merged-vcfs/<outname>_biallelic_truvari_filtered.vcf.gz ``

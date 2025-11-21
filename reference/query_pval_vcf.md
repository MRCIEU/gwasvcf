# Query based on p-value threshold from vcf

Query based on p-value threshold from vcf

## Usage

``` r
query_pval_vcf(pval, vcf, id = NULL)
```

## Arguments

- pval:

  P-value threshold (NOT -log10)

- vcf:

  VCF object (e.g. from readVcf)

- id:

  If multiple GWAS datasets in the vcf file, the name (sample ID) from
  which to perform the filter

## Value

VCF object

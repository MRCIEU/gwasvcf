# Query pval from vcf file

Query pval from vcf file

## Usage

``` r
query_pval_file(pval, vcffile, id = NULL, build = "GRCh37")
```

## Arguments

- pval:

  P-value threshold (NOT -log10)

- vcffile:

  Path to tabix indexed vcf file

- id:

  If multiple GWAS datasets in the vcf file, the name (sample ID) from
  which to perform the filter

- build:

  Default="GRCh37"

## Value

VCF object

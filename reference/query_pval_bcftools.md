# Query p-value using bcftools

Query p-value using bcftools

## Usage

``` r
query_pval_bcftools(pval, vcffile, id = NULL)
```

## Arguments

- pval:

  P-value threshold (NOT -log10)

- vcffile:

  Path to .vcf.gz GWAS summary data file

- id:

  If multiple GWAS datasets in the vcf file, the name (sample ID) from
  which to perform the filter

## Value

vcf object

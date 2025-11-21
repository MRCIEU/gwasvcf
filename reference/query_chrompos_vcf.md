# Query chrompos from vcf object

Query chrompos from vcf object

## Usage

``` r
query_chrompos_vcf(chrompos, vcf, id = NULL)
```

## Arguments

- chrompos:

  Either vector of chromosome and position ranges e.g. "1:1000" or
  "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.

- vcf:

  VCF object (e.g. from readVcf)

- id:

  If multiple GWAS datasets in the vcf file, the name (sample ID) from
  which to perform the filter

## Value

VCF object

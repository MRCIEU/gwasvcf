# Query chromosome and position using bcftools

Query chromosome and position using bcftools

## Usage

``` r
query_chrompos_bcftools(chrompos, vcffile, id = NULL)
```

## Arguments

- chrompos:

  Either vector of chromosome and position ranges e.g. "1:1000" or
  "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.

- vcffile:

  Path to .vcf.gz GWAS summary data file

- id:

  If multiple GWAS datasets in the vcf file, the name (sample ID) from
  which to perform the filter

## Value

vcf object

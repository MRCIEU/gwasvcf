# Query vcf file, extracting by chromosome and position

Query vcf file, extracting by chromosome and position

## Usage

``` r
query_chrompos_file(chrompos, vcffile, id = NULL, build = "GRCh37")
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

- build:

  Default="GRCh37" Build of vcffile

## Value

VCF object

# Query vcf file, extracting by rsid

Query vcf file, extracting by rsid

## Usage

``` r
query_rsid_file(rsid, vcffile, id = NULL, build = "GRCh37")
```

## Arguments

- rsid:

  Vector of rsids. Use DBSNP build (???)

- vcffile:

  Path to .vcf.gz GWAS summary data file

- id:

  If multiple GWAS datasets in the vcf file, the name (sample ID) from
  which to perform the filter

- build:

  Default="GRCh37" Build of vcffile

## Value

VCF object

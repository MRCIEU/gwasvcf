# Query

Query

## Usage

``` r
query_rsid_bcftools(rsid, vcffile, id = NULL)
```

## Arguments

- rsid:

  Vector of rsids

- vcffile:

  Path to .vcf.gz GWAS summary data file

- id:

  If multiple GWAS datasets in the vcf file, the name (sample ID) from
  which to perform the filter

## Value

VCF object

# Query rsid from vcf object

Query rsid from vcf object

## Usage

``` r
query_rsid_vcf(rsid, vcf, id = NULL)
```

## Arguments

- rsid:

  Vector of rsids

- vcf:

  VCF object (e.g. from readVcf)

- id:

  If multiple GWAS datasets in the vcf file, the name (sample ID) from
  which to perform the filter

## Value

VCF object

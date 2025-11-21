# Query rsid from file using rsidx index

See create_rsidx_index

## Usage

``` r
query_rsid_rsidx(rsid, vcffile, id = NULL, rsidx)
```

## Arguments

- rsid:

  Vector of rsids

- vcffile:

  Path to .vcf.gz GWAS summary data file

- id:

  If multiple GWAS datasets in the vcf file, the name (sample ID) from
  which to perform the filter

- rsidx:

  Path to rsidx index file

## Value

vcf object

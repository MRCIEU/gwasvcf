# Query pval from file using pvali index

See create_pvali_index

See create_pvali_index

## Usage

``` r
query_pval_sqlite3(pval, vcffile, id = NULL, pvali)

query_pval_sqlite3(pval, vcffile, id = NULL, pvali)
```

## Arguments

- pval:

  pval threshold

- vcffile:

  Path to .vcf.gz GWAS summary data file

- id:

  If multiple GWAS datasets in the vcf file, the name (sample ID) from
  which to perform the filter

- pvali:

  Path to pval index file

## Value

vcf object

vcf object

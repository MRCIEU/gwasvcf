# Create LD reference sqlite database for tags

This is used for looking up proxies

## Usage

``` r
create_ldref_sqlite(bfile, dbname, tag_r2 = 0.6)
```

## Arguments

- bfile:

  path to plink file

- dbname:

  dbname to produce (overwrites existing if exists)

- tag_r2:

  minimum tag r2

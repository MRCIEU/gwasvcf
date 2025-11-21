# Create new index from existing index using a subset of rsids

Note this requires a modified version of plink that allows ld-window-r2
flag for –r option. Available here:
https://github.com/explodecomputer/plink-ng

## Usage

``` r
create_rsidx_sub_index(rsid, rsidx, newindex)
```

## Arguments

- rsid:

  Vector of rsids

- rsidx:

  Existing index

- newindex:

  New index (Note: will delete existing file if exists)

## Value

NULL, creates new index file

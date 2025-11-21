# Find LD proxies for a set of SNPs

Find LD proxies for a set of SNPs

## Usage

``` r
get_ld_proxies(
  rsid,
  bfile,
  searchspace = NULL,
  tag_kb = 5000,
  tag_nsnp = 5000,
  tag_r2 = 0.6,
  threads = 1,
  out = tempfile()
)
```

## Arguments

- rsid:

  list of rs IDs

- bfile:

  ld reference panel

- searchspace:

  Optional list of rs IDs to use as potential proxies

- tag_kb:

  =5000 Proxy parameter

- tag_nsnp:

  =5000 Proxy parameter

- tag_r2:

  =0.6 Proxy parameter

- threads:

  Number of threads to use (=1)

- out:

  temporary output file

## Value

data frame

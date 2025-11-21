# Extract SNPs from vcf file

Finds proxies if necessary

## Usage

``` r
proxy_match(
  vcf,
  rsid,
  bfile = NULL,
  proxies = "yes",
  tag_kb = 5000,
  tag_nsnp = 5000,
  tag_r2 = 0.6,
  threads = 1,
  rsidx = NULL,
  dbfile = NULL
)
```

## Arguments

- vcf:

  vcf file name

- rsid:

  list of rs IDs

- bfile:

  ld reference panel (plink)

- proxies:

  ="yes" If SNPs are absent then look for proxies (yes) or not (no). Can
  also mask all target SNPs and only return proxies (only), for testing
  purposes

- tag_kb:

  =5000 Proxy parameter

- tag_nsnp:

  =5000 Proxy parameter

- tag_r2:

  =0.6 Proxy parameter

- threads:

  Number of threads to use (=1)

- rsidx:

  Path to rsidx index

- dbfile:

  ld tag database (sqlite)

## Value

data frame

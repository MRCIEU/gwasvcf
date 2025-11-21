# Query data from vcf file

Read in GWAS summary data with filters on datasets (if multiple datasets
per file) and/or chromosome/position, rsids or pvalues. Chooses most
optimal choice for the detected operating system. Typically chrompos
searches are the fastest. On Windows, rsid or pvalue filters from a file
will be slow.

## Usage

``` r
query_gwas(
  vcf,
  chrompos = NULL,
  rsid = NULL,
  pval = NULL,
  id = NULL,
  rsidx = NULL,
  pvali = NULL,
  build = "GRCh37",
  os = Sys.info()[["sysname"]],
  proxies = "no",
  bfile = NULL,
  dbfile = NULL,
  tag_kb = 5000,
  tag_nsnp = 5000,
  tag_r2 = 0.6,
  threads = 1
)
```

## Arguments

- vcf:

  Path or URL to GWAS-VCF file or VCF object e.g. output from
  [`VariantAnnotation::readVcf()`](https://rdrr.io/pkg/VariantAnnotation/man/readVcf-methods.html)
  or
  [`create_vcf()`](https://mrcieu.github.io/gwasvcf/reference/create_vcf.md)

- chrompos:

  Either vector of chromosome and position ranges e.g. "1:1000" or
  "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.

- rsid:

  Vector of rsids

- pval:

  P-value threshold (NOT -log10)

- id:

  If multiple GWAS datasets in the vcf file, the name (sample ID) from
  which to perform the filter

- rsidx:

  Path to rsidx index file

- pvali:

  Path to pval index file

- build:

  ="GRCh37" Build of vcffile

- os:

  The operating system. Default is as detected. Determines the method
  used to perform query

- proxies:

  ="no" If SNPs are absent then look for proxies (yes) or not (no). Can
  also mask all target SNPs and only return proxies (only), for testing
  purposes. Currently only possible if querying rsid.

- bfile:

  =path to plink bed/bim/fam ld reference panel

- dbfile:

  =path to sqlite tag snp database

- tag_kb:

  =5000 Proxy parameter

- tag_nsnp:

  =5000 Proxy parameter

- tag_r2:

  =0.6 Proxy parameter

- threads:

  =1 NUmber of threads

## Value

vcf object

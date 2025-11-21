# Create pval index from GWAS-VCF file

Create a separate file called `<id>.pvali` which is used to speed up
p-value queries.

## Usage

``` r
create_pval_index_from_vcf(vcffile, maximum_pval, indexname)
```

## Arguments

- vcffile:

  VCF filename

- maximum_pval:

  Maximum p-value to include. Default = 0.05

- indexname:

  index file name to create. Deletes existing file if exists.

# Convert vcf format to granges format

Convert vcf format to granges format

## Usage

``` r
vcf_to_granges(vcf, id = NULL)
```

## Arguments

- vcf:

  Output from readVcf

- id:

  Only accepts one ID, so specify here if there are multiple GWAS
  datasets in the vcf

## Value

GRanges object

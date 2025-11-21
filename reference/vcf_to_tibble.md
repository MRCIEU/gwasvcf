# Convert vcf format to tibble (data frame)

Convert vcf format to tibble (data frame)

## Usage

``` r
vcf_to_tibble(vcf, id = NULL)
```

## Arguments

- vcf:

  Output from readVcf

- id:

  Only accepts one ID, so specify here if there are multiple GWAS
  datasets in the vcf

## Value

GRanges object

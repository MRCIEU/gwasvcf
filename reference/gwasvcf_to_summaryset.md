# Create a SummarySet

Returns a gwasglue2 SummarySet object

## Usage

``` r
gwasvcf_to_summaryset(vcf)
```

## Arguments

- vcf:

  Path or URL to GWAS-VCF file or VCF object e.g. output from
  [`VariantAnnotation::readVcf()`](https://rdrr.io/pkg/VariantAnnotation/man/readVcf-methods.html),
  [`create_vcf()`](https://mrcieu.github.io/gwasvcf/reference/create_vcf.md)
  or
  [`query_gwas()`](https://mrcieu.github.io/gwasvcf/reference/query_gwas.md)

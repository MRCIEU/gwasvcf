# Changelog

## gwasvcf 0.1.5

- [`proxy_match()`](https://mrcieu.github.io/gwasvcf/reference/proxy_match.md)
  is now robust to VCFs with multiple samples (thanks
  [@jwr-git](https://github.com/jwr-git))

## gwasvcf 0.1.4

- Add sqlite3 to DESCRIPTION SystemRequirements for
  [`create_pval_index_from_vcf()`](https://mrcieu.github.io/gwasvcf/reference/create_pval_index_from_vcf.md)
- Update some URLs within the package documentation

## gwasvcf 0.1.3

- Fix for security message in
  [`get_ld_proxies()`](https://mrcieu.github.io/gwasvcf/reference/get_ld_proxies.md)
  (thanks [@mattlee821](https://github.com/mattlee821))

## gwasvcf 0.1.2

- New
  [`gwasvcf_to_summaryset()`](https://mrcieu.github.io/gwasvcf/reference/gwasvcf_to_summaryset.md)
  function to create a [gwasglue2](https://mrcieu.github.io/gwasglue2/)
  SummarySet object from a vcf file
- Fixed error in
  [`get_ld_proxies()`](https://mrcieu.github.io/gwasvcf/reference/get_ld_proxies.md)
  related with argument `validate`, deprecated in
  [`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
  (tibble 2.0.0)

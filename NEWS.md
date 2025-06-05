# gwasvcf 0.1.4

* Add sqlite3 to DESCRIPTION SystemRequirements for `create_pval_index_from_vcf()`

# gwasvcf 0.1.3

* Fix for security message in `get_ld_proxies()` (thanks @mattlee821)

# gwasvcf 0.1.2

* New `gwasvcf_to_summaryset()` function to create a [gwasglue2](https://mrcieu.github.io/gwasglue2/) SummarySet object from a vcf file
* Fixed error in `get_ld_proxies()` related with argument `validate`, deprecated in `as_tibble()` (tibble 2.0.0)

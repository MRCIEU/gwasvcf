context("Querying vcf files")
library(gwasvcftools)

test_that("parse_chrompos", {
	expect_equal(parse_chrompos("1:10000") %>% length, 1)
	expect_equal(parse_chrompos("1:10000-100000") %>% length, 1)
	expect_equal(parse_chrompos(c("1:10000-10000", "2:100-200")) %>% length, 2)
	expect_equal(parse_chrompos(tibble(chrom=c(1,2),start=c(10000,100), end=c(10000,200))) %>% length, 2)
})




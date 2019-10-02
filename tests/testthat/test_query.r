context("Querying vcf files")
library(gwasvcftools)

fn <- system.file("data","data.vcf.gz", package="gwasvcftools")
vcf <- readVcf(fn)


test_that("parse_chrompos", {
	expect_equal(parse_chrompos("1:10000") %>% length, 1)
	expect_equal(parse_chrompos("1:10000-100000") %>% length, 1)
	expect_equal(parse_chrompos(c("1:10000-10000", "2:100-200")) %>% length, 2)
	expect_equal(parse_chrompos(tibble(chrom=c(1,2),start=c(10000,100), end=c(10000,200))) %>% length, 2)
})


test_that("granges_from_vcf", {
	g <- granges_from_vcf(vcf)
	expect_equal(length(g), length(vcf))
})


test_that("query_chrompos_file", {
	g <- parse_chrompos("1:800000-1000000")
	v <- query_chrompos_file(g, fn)
	expect_equal(length(v), 3)
})


test_that("query_rsid_file", {
	v <- query_rsid_file(c("rs3128126", "rs3121561", "rs3813193"), fn)
	expect_equal(length(v), 3)
})


test_that("query_pval_file", {
	v <- query_pval_file(0.2, fn)
	expect_true(length(v) < 92)
	expect_true(length(v) > 5)
})


test_that("query_chrompos_vcf", {
	v <- query_chrompos_vcf("1:800000-1000000", vcf)
	expect_equal(length(v), 3)
})


test_that("query_rsid_vcf", {
	v <- query_rsid_vcf(c("rs3128126", "rs3121561", "rs3813193"), vcf)
	expect_equal(length(v), 3)
})


test_that("query_pval_vcf", {
	v <- query_pval_vcf(0.2, vcf)
	expect_true(length(v) < 92)
	expect_true(length(v) > 5)
})


test_that("query_rsid_bcftools", {
	v <- query_rsid_bcftools(c("rs3128126", "rs3121561", "rs3813193"), fn)
	expect_equal(length(v), 3)
})


test_that("query_pval_bcftools", {
	v <- query_pval_bcftools(0.2, fn)
	expect_true(length(v) < 92)
	expect_true(length(v) > 5)
})


test_that("query_chrompos_vcf", {
	v <- query_chrompos_bcftools("1:800000-1000000", fn)
	expect_equal(length(v), 3)
})


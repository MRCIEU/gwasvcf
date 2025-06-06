context("Querying vcf files with pval index")
library(gwasvcf)

fn <- system.file("extdata","data.vcf.gz", package="gwasvcf")
vcf <- VariantAnnotation::readVcf(fn)
if (Sys.info()["sysname"] != "Windows") set_bcftools()

indexname <- tempfile()

test_that("create index", {
  skip_on_os(c("windows", "linux"))
	create_pval_index_from_vcf(fn, 0.4, indexname)
	expect_true(file.exists(indexname))
})

test_that("read in", {
  skip_on_os(c("windows", "linux"))
	out <- query_pvali(0.05, indexname)
	expect_equal(nrow(out), 7)
})

test_that("query with pvali", {
  skip_on_os(c("windows", "linux"))
	b <- query_gwas(fn, pval=0.05, pvali=indexname)
	expect_equal(nrow(b), 7)
})

test_that("query with pvali", {
  skip_on_os("windows")
	b <- query_gwas(fn, pval=0.05)
	expect_equal(nrow(b), 7)
})

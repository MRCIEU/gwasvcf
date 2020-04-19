context("Querying vcf files with rsidx")
library(gwasvcf)


fn <- system.file("extdata","data.vcf.gz", package="gwasvcf")
vcf <- VariantAnnotation::readVcf(fn)

indexname <- tempfile()

test_that("create index", {
	create_rsidx_index_from_vcf(fn, indexname)
	expect_true(file.exists(indexname))
})

test_that("read in", {
	out <- query_rsidx(head(names(vcf)), indexname)
	expect_true(nrow(out) == 6)
})

test_that("create sub index", {
	newname <- tempfile()
	create_rsidx_sub_index(head(names(vcf)), indexname, newname)
	expect_true(file.exists(newname))
})

test_that("query with rsidx", {
	a <- query_gwas(fn, rsid=head(names(vcf)))
	b <- query_gwas(fn, rsid=head(names(vcf)), rsidx=indexname)
	expect_true(all(names(a) == names(b)))
})


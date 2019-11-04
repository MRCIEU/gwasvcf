context("coloc")
library(gwasvcf)

fn <- system.file("extdata","data.vcf.gz", package="gwasvcf")
vcf1 <- VariantAnnotation::readVcf(fn)
vcf2 <- VariantAnnotation::readVcf(fn)

test_that("coloc",
	a <- perform_coloc(vcf1, vcf2)
	expect_true(is.list(a))
)


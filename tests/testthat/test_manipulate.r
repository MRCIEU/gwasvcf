context("VCF manipulations")
library(gwasvcf)

fn <- system.file("extdata","data.vcf.gz", package="gwasvcf")
vcf1 <- VariantAnnotation::readVcf(fn)[1:70,]
vcf2 <- VariantAnnotation::readVcf(fn)[40:90,]
vcf3 <- VariantAnnotation::readVcf(fn)[60:92,]
vcf4 <- VariantAnnotation::readVcf(fn)[65:92,]
l <- list(vcf1, vcf2, vcf3, vcf4)
set_bcftools()

# Need to check what happens with multiallelic variants

test_that("vcflist_overlaps", {
	o <- vcflist_overlaps(vcflist=list(vcf1, vcf2), chrompos=NULL)
	expect_true(all(sapply(o, length) == 31) & length(o) == 2)

	o <- vcflist_overlaps(vcflist=list(vcf1, vcf2, vcf3, fn), chrompos="1:1-10000000")
	expect_true(all(sapply(o, length) == 11) & length(o) == 4)

	o <- vcflist_overlaps(vcflist=list(vcf1, vcf2, vcf3, vcf4), chrompos="1:1-10000000")
	expect_true(all(sapply(o, length) == 6) & length(o) == 4)

	o <- vcflist_overlaps(vcflist=list(fn, fn), chrompos="1:1-10000000")
	expect_true(all(sapply(o, length) == 92) & length(o) == 2)

	o <- vcflist_overlaps(vcflist=list(fn, fn), chrompos="2:1-10000000")
	expect_true(all(sapply(o, length) == 0) & length(o) == 2)	
})





context("Getting LD proxies")
library(gwasvcf)
library(genetics.binaRies)

vcffile <- system.file("extdata","data.vcf.gz", package="gwasvcf")
vcf <- VariantAnnotation::readVcf(vcffile)
bfile <- system.file("extdata","eur.bed", package="gwasvcf") %>% gsub(".bed", "", .)

set_plink()

test_that("query native", {

	set_bcftools(NULL)
	a <- query_gwas(vcffile, rsid="rs4970420")
	expect_equal(nrow(a), 1)

	a <- query_gwas(vcf, rsid="rs4970420")
	expect_equal(nrow(a), 1)

	a <- query_gwas(vcffile, rsid="rs4442317")
	expect_equal(nrow(a), 0)

	a <- query_gwas(vcf, rsid="rs4442317")
	expect_equal(nrow(a), 0)

	a <- query_gwas(vcffile, rsid="rs4442317", proxies="yes", bfile=bfile, tag_r2=0.05)
	expect_equal(nrow(a), 1)

	a <- query_gwas(vcf, rsid="rs4442317", proxies="yes", bfile=bfile, tag_r2=0.05)
	expect_equal(nrow(a), 1)

	a <- query_gwas(vcffile, rsid="rs9729550", proxies="only", bfile=bfile, tag_r2=0.05)
	expect_equal(nrow(a), 1)

	a <- query_gwas(vcf, rsid="rs9729550", proxies="only", bfile=bfile, tag_r2=0.05)
	expect_equal(nrow(a), 1)

})


test_that("query bcftools", {

	set_bcftools()
	a <- query_gwas(vcffile, rsid="rs4970420")
	expect_equal(nrow(a), 1)

	a <- query_gwas(vcf, rsid="rs4970420")
	expect_equal(nrow(a), 1)

	a <- query_gwas(vcffile, rsid="rs4442317")
	expect_equal(nrow(a), 0)

	a <- query_gwas(vcf, rsid="rs4442317")
	expect_equal(nrow(a), 0)

	a <- query_gwas(vcffile, rsid="rs4442317", proxies="yes", bfile=bfile, tag_r2=0.05)
	expect_equal(nrow(a), 1)

	a <- query_gwas(vcffile, rsid=c("rs12565286","rs4442317"), proxies="yes", bfile=bfile, tag_r2=0.05)
	expect_equal(nrow(a), 2)

	a <- query_gwas(vcf, rsid="rs4442317", proxies="yes", bfile=bfile, tag_r2=0.05)
	expect_equal(nrow(a), 1)

	a <- query_gwas(vcf, rsid=c("rs12565286","rs4442317"), proxies="yes", bfile=bfile, tag_r2=0.05)
	expect_equal(nrow(a), 2)

	a <- query_gwas(vcffile, rsid="rs9729550", proxies="only", bfile=bfile, tag_r2=0.05)
	expect_equal(nrow(a), 1)

	a <- query_gwas(vcf, rsid="rs9729550", proxies="only", bfile=bfile, tag_r2=0.05)
	expect_equal(nrow(a), 1)

})

test_that("alignment native", {
	set_bcftools(NULL)
	rsid <- names(SummarizedExperiment::rowRanges(vcf))
	a <- proxy_match(vcf, rsid, bfile, proxies="only")
	b <- query_gwas(vcf, rsid=rsid)
	index <- match(names(b), names(a))
	names(b) == names(a)[index]
	expect_true(cor(vcf_to_granges(b)$ES, vcf_to_granges(a)$ES[index], use="pair") > 0.5)
})

test_that("alignment bcftools", {
	set_bcftools()
	rsid <- names(SummarizedExperiment::rowRanges(vcf))
	a <- proxy_match(vcf, rsid, bfile, proxies="only")
	b <- query_gwas(vcf, rsid=rsid)
	index <- match(names(b), names(a))
	names(b) == names(a)[index]
	expect_true(cor(vcf_to_granges(b)$ES, vcf_to_granges(a)$ES[index], use="pair") > 0.5)
})

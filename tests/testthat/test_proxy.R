context("Getting LD proxies")
library(gwasvcftools)

vcffile <- system.file("data","data.vcf.gz", package="gwasvcftools")
vcf <- readVcf(vcffile)
bfile <- system.file("data","eur.bed", package="gwasvcftools") %>% gsub(".bed", "", .)

test_that("query", {
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

test_that("alignment", {
	rsid <- names(rowRanges(vcf))
	a <- proxy_match(vcf, rsid, bfile, proxies="only")
	b <- query_gwas(vcf, rsid=rsid)
	index <- match(names(b), names(a))
	names(b) == names(a)[index]
	expect_true(cor(vcf_to_granges(b)$ES, vcf_to_granges(a)$ES[index], use="pair") > 0.5)
})

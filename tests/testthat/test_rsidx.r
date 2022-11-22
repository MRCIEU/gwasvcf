context("Querying vcf files with rsidx")
library(gwasvcf)


fn <- system.file("extdata","data.vcf.gz", package="gwasvcf")
vcf <- VariantAnnotation::readVcf(fn)

indexname <- tempfile()

test_that("create index", {
  skip_on_os("windows")
	create_rsidx_index_from_vcf(fn, indexname)
	expect_true(file.exists(indexname))
})

test_that("read in", {
  skip_on_os("windows")
	out <- query_rsidx(head(names(vcf)), indexname)
	expect_true(nrow(out) == 6)
})

test_that("create sub index", {
  skip_on_os("windows")
	newname <- tempfile()
	create_rsidx_sub_index(head(names(vcf)), indexname, newname)
	expect_true(file.exists(newname))
})

test_that("query with rsidx", {
  skip_on_os("windows")
	a <- query_gwas(fn, rsid=head(names(vcf)))
	b <- query_gwas(fn, rsid=head(names(vcf)), rsidx=indexname)
	expect_true(all(names(a) == names(b)))
})



fn <- system.file("extdata", "eur.bed", package="gwasvcf") %>% gsub("eur.bed", "eur", .)
dbfile <- tempfile()

set_plink()

test_that("tag db", {
  skip_on_os("windows")
	create_ldref_sqlite(fn, dbfile, 0.04)
	expect_true(file.exists(dbfile))
})

test_that("sqlite_ld_proxies", {
  skip("TODO: check this test")
	m <- data.table::fread(paste0(fn, ".bim")) %>% {sample(.$V2, 100, replace=FALSE)}
	ld <- sqlite_ld_proxies(m, dbfile, 0.2)
	# Requires an expect_* condition here
})

test_that("sqlite proxy", {
  skip("TODO: check this test")
	vcffile <- system.file("extdata","data.vcf.gz", package="gwasvcf")
	set_bcftools()
	a <- query_gwas(vcffile, rsid="rs4442317", proxies="yes", dbfile=dbfile, tag_r2=0.05)
	expect_equal(nrow(a), 1)
})

test_that("sqlite proxy", {
  skip("TODO: check this test")
	vcffile <- system.file("extdata","data.vcf.gz", package="gwasvcf")
	set_bcftools()
	a <- query_gwas(vcffile, rsid=c("rs12565286","rs4442317"), proxies="yes", dbfile=dbfile, tag_r2=0.05)
	expect_equal(nrow(a), 2)
})

test_that("sqlite proxy only", {
  skip("TODO: check this test")
	vcffile <- system.file("extdata","data.vcf.gz", package="gwasvcf")
	set_bcftools()
	a <- query_gwas(vcffile, rsid=c("rs12565286","rs4442317"), proxies="only", dbfile=dbfile, tag_r2=0.05)
	expect_equal(nrow(a), 2)
	b <- a %>% vcf_to_tibble
	expect_true(all(! b$ID %in% c("rs12565286","rs4442317")))
})

test_that("sqlite proxy no result", {
  skip_on_os("windows")
  vcffile <- system.file("extdata","data.vcf.gz", package="gwasvcf")
	set_bcftools()
	a <- query_gwas(vcffile, rsid="rs4442317", proxies="yes", dbfile=dbfile, tag_r2=0.5)
	expect_equal(nrow(a), 0)
})

unlink(dbfile)

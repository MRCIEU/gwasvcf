context("Creating and harmonising data")
library(gwasvcf)
library(magrittr)

fn <- system.file("data","data.vcf.gz", package="gwasvcf")
V <- readVcf(fn) 
vv <- V %>% vcf_to_granges %>% as_tibble

test_that("create vcf", {
	out <- vv %$% create_vcf(chrom=seqnames, pos=start, nea=REF, ea=ALT, snp=ID, ea_af=AF, effect=ES, se=SE, pval=10^-LP, n=SS, name="a")
	writeVcf(out, file="temp.vcf")
	expect_true(file.exists("temp.vcf"))
})




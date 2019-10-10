context("Getting LD proxies")
library(gwasvcftools)



vcf <- "~/data/gwas/IEU-a-2.vcf.gz"
bfile <- "~/data/ld_files/data_maf0.01_rs_ref"
rsid <- TwoSampleMR::extract_instruments(2)$SNP

a <- proxy_match(vcf, rsid, bfile, proxies="only")
b <- query_gwas(vcf, rsid=rsid)
index <- match(names(b), names(a))
names(b) == names(a)[index]
all(sign(vcf_to_granges(b)$ES) == sign(vcf_to_granges(a)$ES[index]), na.rm=T)


rsid <- c("rs12659333","rs13328017","rs28616030","rs60199041","rs2963074", "rs61980157","rs4864550","rs1396028","rs6486162","rs56102775", "rs147785482","rs7044202","rs619228","rs5755921", "rs377394")
a <- proxy_match(vcf, rsid, bfile, proxies="yes")
b <- query_gwas(vcf, rsid=rsid)
index <- match(names(b), names(a))
names(b) == names(a)[index]
all(sign(vcf_to_granges(b)$ES) == sign(vcf_to_granges(a)$ES[index]), na.rm=T)



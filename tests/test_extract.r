load_all()
s <- fread("~/mr-eve/mr-eve/instruments.txt")
a <- extract(
	bcf="~/mr-eve/gwas-files/7/data.bcf", 
	snplist=s, 
	tempname="temp", 
	proxies="yes", 
	bfile="~/mr-eve/vcf-reference-datasets/1000g_filtered/data_maf0.01_rs_snps",
	vcf="~/mr-eve/vcf-reference-datasets/1000g/1kg_v3_nomult.bcf"
)

a <- extract(
	bcf="~/mr-eve/gwas-files/7/data.bcf", 
	snplist=s[[6]], 
	tempname="temp", 
	proxies="yes", 
	bfile="~/mr-eve/vcf-reference-datasets/1000g_filtered/data_maf0.01_rs_snps",
	vcf="~/mr-eve/vcf-reference-datasets/1000g/1kg_v3_nomult.bcf"
)


a <- extract(
	bcf="~/mr-eve/gwas-files/2/data.bcf", 
	snplist=s[[6]], 
	tempname="temp", 
	proxies="yes", 
	bfile="~/mr-eve/vcf-reference-datasets/1000g_filtered/data_maf0.01_rs_snps",
	vcf="~/mr-eve/vcf-reference-datasets/1000g/1kg_v3_nomult.bcf"
)

a <- extract(
	bcf="~/mr-eve/gwas-files/2/data.bcf", 
	snplist=s, 
	tempname="temp", 
	proxies="yes", 
	bfile="~/mr-eve/vcf-reference-datasets/1000g_filtered/data_maf0.01_rs_snps",
	vcf="~/mr-eve/vcf-reference-datasets/1000g/1kg_v3_nomult.bcf"
)




library(devtools)
load_all()
a <- TwoSampleMR::extract_instruments(2)
fn <- system.file("data","IEU-a-2.vcf.gz", package="gwasvcftools")
ldref <- "~/repo/mr-base-api/app/ld_files/data_maf0.01_rs"

o <- get_ld_proxies(a$SNP, fn, ldref, tempfile())

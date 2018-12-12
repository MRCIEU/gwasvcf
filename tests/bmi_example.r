if(!require(gwasvcftools))
{
	if(!required(devtools)) install.packages("devtools")
	devtools::install_github("MRCIEU/gwasvcftools")
}
library(gwasvcftools)
library(argparse)

# create parser object
parser <- ArgumentParser()
parser$add_argument('--snplist', required=TRUE)
parser$add_argument('--bcf-dir', required=TRUE)
parser$add_argument('--gwas-id', required=TRUE)
parser$add_argument('--out', required=TRUE)
parser$add_argument('--bfile', required=TRUE)
parser$add_argument('--get-proxies', action="store_true", default=FALSE)
parser$add_argument('--vcf-ref', required=FALSE)
parser$add_argument('--tag-r2', type="double", default=0.6)
parser$add_argument('--tag-kb', type="double", default=5000)
parser$add_argument('--tag-nsnp', type="double", default=5000)
parser$add_argument('--palindrome-freq', type="double", default=0.4)
parser$add_argument('--no-clean', action="store_true", default=FALSE)
parser$add_argument('--rdsf-config', required=FALSE, default='')
parser$add_argument('--instrument-list', required=FALSE)


# args <- parser$parse_args()
setwd("~/mr-eve/gwas-instrument-subsets/scripts")
args <- parser$parse_args(c("--bfile", "../../vcf-reference-datasets/ukb/ukb_ref", "--gwas-id", "2", "--snplist", "temp.snplist", "--no-clean", "--out", "out", "--bcf-dir", "../../gwas-files", "--vcf-ref", "../../vcf-reference-datasets/1000g/1kg_v3_nomult.bcf", "--get-proxies"))
print(args)
tempname <- tempfile(pattern="extract", tmpdir=dirname(args[['out']]))
bcf <- file.path(args[['bcf_dir']], args[['gwas_id']], "harmonised.bcf")
snplist <- scan(args[['snplist']], what=character())


# Test Different proxy options

o1 <- wrapper(bcf, snplist, tempname, "yes", args[["bfile"]], args[["vcf_ref"]])
dim(o1)
o2 <- wrapper(bcf, snplist, tempname, "no", args[["bfile"]])
dim(o2)
o3 <- wrapper(bcf, snplist, tempname, "only", args[["bfile"]], args[["vcf_ref"]])
dim(o3)


# Check that proxies are correctly oriented
# Expect to see that the proxies (o3) have effect sizes that strongly correlate with the true effect sizes (o2)

a <- merge(o3, o2, by="ID")
i <- a$ALT.x == a$ALT.y
table(i)
cor(a$B.x, a$B.y)
plot(a$B.x, a$B.y)


# Finally, check that the original elastic files are on the same strand as the harmonised data

o <- fread("gunzip -c ../../gwas-files/2/elastic.gz", he=FALSE)
temp <- merge(o3, o, by.x="ID", by.y="V1")
dim(temp)
i <- temp$REF != temp$V2
table(i)
cor(temp$B, temp$V5)
temp$B[i] <- temp$B[i] * -1
cor(temp$B, temp$V5)
plot(temp$B, temp$V5)



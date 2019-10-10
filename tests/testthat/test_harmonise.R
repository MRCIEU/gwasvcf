context("Creating and harmonising data")
library(gwasvcftools)

fn <- system.file("data","data.vcf.gz", package="gwasvcftools")
V <- readVcf(fn) 
vv <- V %>% vcf_to_granges %>% as_tibble

out <- vv %$% create_vcf(chrom=seqnames, pos=start, nea=REF, ea=ALT, snp=ID, ea_af=AF, effect=ES, se=SE, pval=10^-LP, n=SS, name="a")
writeVcf(out, file="temp.vcf")



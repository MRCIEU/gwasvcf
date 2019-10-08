#' Extract SNPs from vcf 
#'
#'
#' @param snplist list of rs IDs or table of chrom and pos
#' @param vcf vcf file name
#' @param out temporary filename
#' @param sel='\%CHROM \%POS \%ID \%REF \%ALT \%EFFECT \%SE \%L10PVAL \%N \%AF\n' What to extract from vcf
#' @param logpval Whether to return -log10 or standard p-value. Default is standard (though it is stored as logged in vcf to retain precision)
#'
#' @export
#' @return data frame
extract_from_vcf <- function(snplist, vcf, out, sel='%CHROM %POS %ID %REF %ALT [%ES %SE %LP %SS %AF]\n', logpval=FALSE)
{
	require(data.table)
	snplistname <- paste0(out, ".snplist")
	extractname <- paste0(out, ".extract")

	if(is.vector(snplist))
	{
		write.table(snplist, file=snplistname, row=F, col=F, qu=F, sep="\t")
		nsnp <- length(snplist)
		cmd <- stringf("%s view -i'ID=@%s.snplist' %s | %s query -f'%s' > %s", get_bcftools_binary(), out, vcf, get_bcftools_binary(), sel, extractname)
	} else {
		stopifnot(is.data.frame(snplist))
		write.table(snplist[,1:2], file=snplistname, row=F, col=F, qu=F, sep="\t")
		nsnp <- nrow(snplist)
		cmd <- stringf("%s query -R %s -f'%s' %s > %s", get_bcftools_binary(), snplistname, sel, vcf, extractname)
	}
	system(cmd)
	o <- data.table::fread(extractname, header=FALSE, sep=" ", na.strings=".")
	if(logpval)
	{
		which({gsub("[[:space:]]", "", sel) %>% strsplit(sel, split="%")}[[1]] == "LP") - 1
		o[,8] <- 10^-o[,8]
	}
	message("Extracted ", nrow(o), " out of ", nsnp, " SNPs")
	unlink(extractname)
	unlink(snplistname)
	return(o)
}


#' Find LD proxies for a set of SNPs
#'
#' @param rsid list of rs IDs
#' @param vcf vcf file name
#' @param bfile ld reference panel
#' @param out temporary output file
#' @param tag_kb=5000 Proxy parameter
#' @param tag_nsnp=5000 Proxy parameter
#' @param tag_r2=0.6 Proxy parameter
#'
#' @export
#' @return data frame
get_ld_proxies <- function(rsid, vcf, bfile, out=tempfile(), tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1)
{
	require(dplyr)
	require(data.table)
	searchspacename <- paste0(out, ".searchspace")
	searchspacename1 <- paste0(out, ".searchspace1")
	targetsname <- paste0(out, ".targets")
	outname <- paste0(out, ".targets.ld.gz")

	cmd <- paste0("bcftools query -f'%ID\n' ", vcf, " > ", searchspacename)
	system(cmd)
	write.table(rsid, file=targetsname, row=FALSE, col=FALSE, qu=FALSE)
	cmd <- paste0("cat ", targetsname, " ", searchspacename, " > ", searchspacename1)
	system(cmd)

	cmd <- paste0(get_plink_binary(),
		" --bfile ", bfile, 
		" --extract ", searchspacename1,
		" --keep-allele-order ",
		" --r in-phase with-freqs gz",
		" --ld-snp-list ", targetsname,
		" --ld-window-kb ", tag_kb,
		" --ld-window ", tag_nsnp,
		" --out ", targetsname,
		" --threads ", threads
	)
	system(cmd)

	ld <- data.table::fread(paste0("gunzip -c ", outname), header=TRUE) %>%
		dplyr::filter(R^2 > tag_r2) %>%
		dplyr::filter(SNP_A != SNP_B) %>%
		dplyr::mutate(PHASE=gsub("/", "", PHASE))
	temp <- do.call(rbind, strsplit(ld$PHASE, "")) %>% as_tibble
	names(temp) <- c("A1", "B1", "A2", "B2")
	ld <- cbind(ld, temp) %>% as_tibble()
	ld <- arrange(ld, desc(abs(R))) %>%
		dplyr::filter(!duplicated(SNP_A))
	unlink(searchspacename)
	unlink(searchspacename1)
	unlink(targetsname)
	unlink(paste0(targetsname, c(".log", ".nosex")))
	unlink(outname)
	return(ld)
	index <- ld$R < 0
	temp <- ld$B1[index]
	ld$B1[index] <- ld$B2[index]
	ld$B2[index] <- temp
	ld$R[index] <- ld$R[index] * -1
	return(ld)
}


#' Extract SNPs from vcf file
#'
#' Finds proxies if necessary
#'
#' @param vcf vcf file name
#' @param rsid list of rs IDs
#' @param bfile ld reference panel
#' @param proxies="yes" If SNPs are absent then look for proxies (yes) or not (no). Can also mask all target SNPs and only return proxies (only), for testing purposes
#' @param vcf_ref Optional, if provided then final dataset is aligned to reference
#' @param tag_kb=5000 Proxy parameter
#' @param tag_nsnp=5000 Proxy parameter
#' @param tag_r2=0.6 Proxy parameter
#' @param threads=1
#'
#' @export
#' @return data frame
proxy_match <- function(vcf, rsid, bfile, proxies="yes", tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1)
{

	if(proxies=="yes")
	{
		o <- query_gwas(vcf, rsid=rsid)
		missing <- rsid[!rsid %in% names(o)]
		if(length(missing) != 0)
		{
			ld <- get_ld_proxies(missing, vcf, bfile, tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1)
		} else {
			return(o)
		}
	} else if(proxies == "only") {
		ld <- get_ld_proxies(rsid, vcf, bfile, tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1)
	} else if(proxies == "no") {
		o <- query_gwas(vcf, rsid=rsid)
		return(o)
	} else {
		stop('proxies must be "yes", "no" or "only"')
	}

	e <- query_gwas(vcf, rsid=ld$SNP_B)
	e <- e[names(e) %in% ld$SNP_B, ]
	index <- match(names(e), ld$SNP_B)
	ld <- ld[index,]
	stopifnot(all(ld$SNP_B == names(e)))
	sign_index <- rowRanges(e)$REF == ld$B1

	gr <- GRanges(ld$CHR_A, IRanges(start=ld$BP_A, end=ld$BP_A, names=ld$SNP_A)) %>% sort
	prox <- VCF(
		rowRanges = gr,
		colData = colData(e),
		info = info(e),
		exptData = list(
			header = header(e), 
			fixed = DataFrame(paramRangeID=as.factor(rep(NA, nrow(ld))),  REF=DNAStringSet(ld$A1), ALT=DNAStringSetList(as.list(ld$A2)), QUAL=as.numeric(NA), FILTER="PASS")
		),
		geno = SimpleList(
			lapply(geno(e), `dimnames<-`, NULL)
		)
	)

	geno(prox)$ES[!sign_index] <- {unlist(geno(prox)$ES[!sign_index]) * -1} %>% as.list
	geno(prox)$PR <- matrix(ld$SNP_B, length(ld$SNP_B), 1)
	geno(header(prox)) <- rbind(geno(header(prox)), 
		DataFrame(row.names="PR", Number="1", Type="String", Description="Proxy rsid")
	)

	if(proxies == "only")
	{
		return(prox)
	} else {
		geno(o)$PR <- matrix(rep(NA, length(o)), length(o), 1)
		geno(header(o)) <- rbind(geno(header(o)), DataFrame(row.names="PR", Number="1", Type="String", Description="Proxy rsid"))
		return(BiocGenerics::rbind(o, prox))
	}
}


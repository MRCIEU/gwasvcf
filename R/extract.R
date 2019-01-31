#' Extract SNPs from vcf 
#'
#'
#' @param snplist list of rs IDs or table of chrom and pos
#' @param bcf bcf file name
#' @param out temporary filename
#' @param sel='%CHROM %POS %ID %REF %ALT %B %SE %L10PVAL %N %AF\n' What to extract from bcf
#' @param logpval Whether to return -log10 or standard p-value. Default is standard (though it is stored as logged in bcf to retain precision)
#'
#' @export
#' @return data frame
extract_from_bcf <- function(snplist, bcf, out, sel='%CHROM %POS %ID %REF %ALT %B %SE %L10PVAL %N %AF\n', logpval=FALSE)
{
	require(data.table)
	snplistname <- paste0(out, ".snplist")
	extractname <- paste0(out, ".extract")

	
	if(is.vector(snplist))
	{
		write.table(snplist, file=snplistname, row=F, col=F, qu=F, sep="\t")
		nsnp <- length(snplist)
		cmd <- paste0("bcftools view -i'ID=@", out, ".snplist' ", bcf, " | bcftools query -f'", sel, "' > ", extractname)
	} else {
		stopifnot(is.data.frame(snplist))
		write.table(snplist[,1:2], file=snplistname, row=F, col=F, qu=F, sep="\t")
		nsnp <- nrow(snplist)
		cmd <- paste0("bcftools query -R ", snplistname, " -f'", sel, "' ", bcf, " > ", extractname)
	}
	system(cmd)
	o <- data.table::fread(extractname, header=FALSE, sep=" ", na.strings=".")
	if(logpval)
	{
		which({gsub("[[:space:]]", "", sel) %>% strsplit(sel, split="%")}[[1]] == "L10PVAL") - 1
		o[,8] <- 10^-o[,8]
	}
	message("Extracted ", nrow(o), " out of ", nsnp, " SNPs")
	unlink(extractname)
	unlink(snplistname)
	return(o)
}


#' Find LD proxies for a set of SNPs
#'
#' @param rsids list of rs IDs
#' @param bcf vcf file name
#' @param bfile ld reference panel
#' @param out temporary output file
#' @param tag_kb=5000 Proxy parameter
#' @param tag_nsnp=5000 Proxy parameter
#' @param tag_r2=0.6 Proxy parameter
#'
#' @export
#' @return data frame
get_ld_proxies <- function(rsids, bcf, bfile, out, tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1)
{
	require(dplyr)
	require(data.table)
	searchspacename <- paste0(out, ".searchspace")
	searchspacename1 <- paste0(out, ".searchspace1")
	targetsname <- paste0(out, ".targets")
	outname <- paste0(out, ".targets.ld.gz")

	cmd <- paste0("bcftools query -f'%ID\n' ", bcf, " > ", searchspacename)
	system(cmd)
	print(head(rsids))
	write.table(rsids, file=targetsname, row=FALSE, col=FALSE, qu=FALSE)
	cmd <- paste0("cat ", targetsname, " ", searchspacename, " > ", searchspacename1)
	system(cmd)

	cmd <- paste0(
		"plink --bfile ", bfile, 
		" --extract ", searchspacename1,
		" --r2 in-phase with-freqs gz",
		" --ld-snp-list ", targetsname,
		" --ld-window-kb ", tag_kb,
		" --ld-window-r2 ", tag_r2,
		" --ld-window ", tag_nsnp,
		" --out ", targetsname,
		" --threads ", threads
	)
	system(cmd)

	ld <- data.table::fread(paste0("gunzip -c ", outname), header=TRUE) %>%
		dplyr::filter(SNP_A != SNP_B) %>%
		dplyr::mutate(PHASE=gsub("/", "", PHASE))
	temp <- do.call(rbind, strsplit(ld$PHASE, "")) %>% as_data_frame
	names(temp) <- c("A1", "B1", "A2", "B2")
	ld <- cbind(ld, temp) %>% as_data_frame()
	ld <- arrange(ld, desc(R2)) %>%
		dplyr::filter(!duplicated(SNP_A))
	unlink(searchspacename)
	unlink(searchspacename1)
	unlink(targetsname)
	unlink(paste0(targetsname, c(".log", ".nosex")))
	unlink(outname)
	return(ld)
}

#' Align proxies
#'
#' Takes output from get_ld_proxies and extract_from_bcf to make sure effect estimates are relative to the correct proxy alleles
#'
#' @param ld output from get_ld_proxies
#' @param e output from extract_from_bcf
#' @param vcf_ref Optional, if provided then final dataset is aligned to reference
#' @param tempfile temporary files for use for extractions
#'
#' @export
#' @return
align_proxies <- function(ld, e, vcf_ref=NULL, tempfile=NULL)
{
	require(dplyr)
	require(magrittr)
	temp <- merge(e, ld, by.x="V3", by.y="SNP_B")
	if(nrow(temp) == 0)

	temp <- subset(temp, (V4 == B1 & V5 == B2) | (V4 == B2 & V5 == B1))
	sw <- temp$V4 == temp$B1
	temp$V6[sw] <- temp$V6[sw] * -1
	r <- temp %$% data_frame(ID = SNP_A, CHROM = CHR_A, POS = BP_A, REF = A1, ALT = A2, B = V6, SE = V7, PVAL = V8, N = V9, AF = MAF_B, proxy.chrom = CHR_B, proxy.pos = BP_B, proxy.id = V3)
	ref <- extract_from_bcf(subset(r, select=c(CHROM, POS)), vcf_ref, tempfile, '%CHROM %POS %ID %REF %ALT\n')

	a1 <- merge(r, ref, by.x="ID", by.y="V3")
	i <- a1$REF == a1$V4
	a1$B[i] <- a1$B[i] * -1
	tt <- a1$REF[i]
	a1$REF[i] <- a1$ALT[i]
	a1$ALT[i] <- tt
	a1 <- subset(a1, select = -c(V1, V2, V4, V5))
	return(a1)
}

#' Write to neo4j csv format
#'
#' @param x data frame
#' @param basename Filename to write to
#' @param header=FALSE Whether to also write a header
#'
#' @export
#' @return NULL
write_out <- function(x, basename, header=FALSE)
{
	g <- gzfile(basename, "w")
	write.table(x, g, row.names=FALSE, col.names=FALSE, na="", sep=",")
	close(g)
	if(header) write.table(x[0,], file=paste0(basename, "_header.csv"), row.names=FALSE, col.names=TRUE, sep=",")
}

is_palindrome <- function(a1, a2)
{
	(a1 == "A" & a2 == "T") |
	(a1 == "T" & a2 == "A") |
	(a1 == "C" & a2 == "G") |
	(a1 == "G" & a2 == "C")
}

#' Extract SNPs from bcf file
#'
#' Finds proxies if necessary
#'
#' @param bcf vcf file name
#' @param snplist list of rs IDs or data frame of chrom and pos
#' @param tempname Temporary file
#' @param proxies="yes" If SNPs are absent then look for proxies (yes) or not (no). Can also mask all target SNPs and only return proxies (only), for testing purposes
#' @param bfile ld reference panel
#' @param vcf_ref Optional, if provided then final dataset is aligned to reference
#' @param tag_kb=5000 Proxy parameter
#' @param tag_nsnp=5000 Proxy parameter
#' @param tag_r2=0.6 Proxy parameter
#'
#' @export
#' @return
extract <- function(bcf, snplist, tempname, proxies="yes", bfile, vcf_ref, tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1)
{
	if(is.vector(snplist))
	{
		rsid <- TRUE
	} else {
		stopifnot(is.data.frame(snplist))
		rsid <- FALSE
	}
	stopifnot(proxies %in% c("yes", "no", "only"))
	if(proxies == "no")
	{
		o <- extract_from_bcf(snplist, bcf, tempname)
		names(o) <- c("CHROM", "POS", "ID", "REF", "ALT", "B", "SE", "PVAL", "N", "AF")
		o$proxy.chrom <- o$proxy.pos <- o$proxy.id <- NA
		return(o)
	}

	if(proxies == "yes")
	{
		o <- extract_from_bcf(snplist, bcf, tempname)

		if(rsid)
		{
			missing_snps <- snplist[!snplist %in% o$V3]
			ld <- get_ld_proxies(missing_snps, bcf, bfile, tempname, threads=threads)
		} else {
			id1 <- paste(snplist[,1], snplist[,2], snplist[,3], snplist[,4])
			missing_snps <- snplist[!id1 %in% paste(o$V1, o$V2, o$V4, o$V5),]
			ld <- get_ld_proxies(missing_snps[,6], bcf, bfile, tempname, threads=threads)
		}
		e <- extract_from_bcf(ld$SNP_B, bcf, tempname)
		a <- align_proxies(ld, e, vcf_ref, tempname)
		names(o) <- c("CHROM", "POS", "ID", "REF", "ALT", "B", "SE", "PVAL", "N", "AF")
		final <- bind_rows(a,o) %>% arrange(CHROM, POS)
		return(final)
	}

	if(proxies == "only")
	{
		if(rsid)
		{
			ld <- get_ld_proxies(snplist, bcf, bfile, tempname, threads=threads)
		} else {
			ld <- get_ld_proxies(snplist[,6], bcf, bfile, tempname, threads=threads)
		}
		e <- extract_from_bcf(ld$SNP_B, bcf, tempname)
		a <- align_proxies(ld, e, vcf_ref, tempname)
		return(a)
	}
}






#' Find binary for bcftools
#'
#' @param 
#'
#' @export
#' @return Path to bcftools
get_bcftools_binary <- function()
{
	switch(Sys.info()[['sysname']],
		Windows = { stop("Sorry, bcftools binary is not available for Windows at the moment. Use the other native functions for querying, or for faster speeds use this package on Mac or Linux")},
		Linux = { system.file("bin", "bcftools_Linux", package="gwasvcftools") },
		Darwin = { system.file("bin", "bcftools_Darwin", package="gwasvcftools") })
}

#' Find binary for plink
#'
#' @param 
#'
#' @export
#' @return Path to plink
get_plink_binary <- function()
{
	switch(Sys.info()[['sysname']],
		Windows = { system.file("bin", "plink.exe", package="gwasvcftools") },
		Linux = { system.file("bin", "plink_Linux", package="gwasvcftools") },
		Darwin = { system.file("bin", "plink_Darwin", package="gwasvcftools") })
}

granges_from_df <- function(df)
{
	GRanges(seqnames=df$chrom, ranges=IRanges(start=df$start, end=df$end))
}


#' Parse chromosome:position
#'
#' Takes data frame or vector of chromosome position ranges and parses to granges object
#'
#' @param chrompos Either vector of chromosome and position ranges e.g. "1:1000" or "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.
#'
#' @export
#' @return GRanges object
parse_chrompos <- function(chrompos)
{

	if("GRanges" %in% class(chrompos))
	{
		return(chrompos)
	} else if(is.data.frame(chrompos)) {
		stopifnot(is.data.frame(chrompos))
		stopifnot(all(c("chrom", "start", "end") %in% names(chrompos)))
		return(granges_from_df(chrompos))
	} else if(!is.character(chrompos)) {
		stop("chrompos must be data frame with columns chrom, start, end, or character vector of <chr:pos> or <chr:start-end>")
	}

	a <- stringr::str_split(chrompos, ":")
	chrom <- sapply(a, function(x) x[1])
	pos <- sapply(a, function(x) x[2])
	i <- grepl("-", pos)
	temp <- stringr::str_split(pos[i], "-")
	pos1 <- pos
	pos2 <- pos
	pos1[i] <- sapply(temp, function(x) {x[1]})
	pos2[i] <- sapply(temp, function(x) {x[2]})
	pos1 <- as.numeric(pos1)
	pos2 <- as.numeric(pos2)
	return(granges_from_df(data.frame(chrom, start=pos1, end=pos2, stringsAsFactors=FALSE)))
}


#' Convert vcf format to granges format
#'
#' @param vcf Output from readVcf
#' @param id Only accepts one ID, so specify here if there are multiple GWAS datasets in the vcf
#'
#' @export
#' @return GRanges object
granges_from_vcf <- function(vcf, id=NULL)
{
	if(is.null(id))
	{
		id <- samples(header(vcf))
	}
	stopifnot(length(id) == 1)
	out <- VariantAnnotation::expand(vcf) %>% 
	geno() %>%
	as.list() %>%
	lapply(., function(x) x[,id, drop=TRUE]) %>%
	bind_cols()
	a <- rowRanges(vcf)
	values(a) <- cbind(values(a), out)
	values(a)$id <- id

	return(a)
}


#' Query vcf file, extracting by chromosome and position
#'
#'
#' @param chrompos Either vector of chromosome and position ranges e.g. "1:1000" or "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.
#' @param vcffile Path to .vcf.gz GWAS summary data file
#' @param build="hg19" Build of vcffile
#'
#' @export
#' @return VCF object
query_chrompos_file <- function(chrompos, vcffile, id=NULL, build="hg19")
{
	chrompos <- parse_chrompos(chrompos)
	tab <- TabixFile(vcffile)
	vcf <- readVcf(tab, build, param=chrompos)
	return(vcf)
}


#' Query vcf file, extracting by rsid
#'
#' @param rsid Vector of rsids. Use DBSNP build (???)
#' @param vcffile Path to .vcf.gz GWAS summary data file
#' @param build="hg19" Build of vcffile
#'
#' @export
#' @return VCF object
query_rsid_file <- function(rsid, vcffile, build="hg19")
{
	message("Note, this is much slower than searching by chromosome/position (e.g. see query_chrompos_file)")
	vcf <- TabixFile(vcffile)
	fil <- function(x)
	{
		grepl(paste(rsid, collapse="|"), x)
	}

	tempfile <- tempfile()
	filterVcf(vcf, build, tempfile, prefilters=FilterRules(list(fil=fil)), verbose=TRUE)
	o <- readVcf(tempfile)
	unlink(tempfile)

	# Grep isn't matching on exact word so do second pass here
	o <- query_rsid_vcf(rsid, o)
	return(o)
}


#' Query pval from vcf file
#'
#' @param pval P-value threshold (NOT -log10)
#' @param vcffile Path to tabix indexed vcf file
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#' @param build="hg19"
#'
#' @export
#' @return VCF object
query_pval_file <- function(pval, vcffile, id=NULL, build="hg19")
{
	if(is.null(id))
	{
		id <- samples(scanVcfHeader(vcffile))
	}
	stopifnot(length(id) == 1)
	message("Filtering p-value based on id ", id)
	message("Note, this is much slower than searching by chromosome/position (e.g. see query_chrompos_file)")
	vcf <- TabixFile(vcffile)
	fil <- function(x)
	{
		geno(x)$LP[,id,drop=TRUE] > -log10(pval)
	}
	tempfile <- tempfile()
	filterVcf(vcf, build, tempfile, filters=FilterRules(list(fil=fil)), verbose=TRUE)
	o <- readVcf(tempfile)
	unlink(tempfile)
	return(o)
}


#' Query chrompos from vcf object
#'
#' @param chrompos Either vector of chromosome and position ranges e.g. "1:1000" or "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.
#' @param vcf VCF object (e.g. from readVcf)
#'
#' @export
#' @return VCF object
query_chrompos_vcf <- function(chrompos, vcf)
{
	chrompos <- parse_chrompos(chrompos)
	i <- IRanges::findOverlaps(SummarizedExperiment::rowRanges(vcf), chrompos) %>% S4Vectors::queryHits() %>% unique %>% sort
	vcf[i,]
}


#' Query rsid from vcf object
#'
#' @param rsid Vector of rsids
#' @param vcf VCF object (e.g. from readVcf)
#'
#' @export
#' @return VCF object
query_rsid_vcf <- function(rsid, vcf)
{
	vcf[rownames(vcf) %in% rsid,]
}


#' Query based on p-value threshold from vcf
#'
#' @param pval P-value threshold (NOT -log10)
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#' @param vcf VCF object (e.g. from readVcf)
#'
#' @export
#' @return VCF object
query_pval_vcf <- function(pval, vcf, id=NULL)
{
	if(is.null(id))
	{
		id <- samples(header(vcf))
	}
	stopifnot(length(id) == 1)
	colid <- which(samples(header(vcf)) == id)
	vcf[geno(vcf)$LP[,colid,drop=TRUE] > -log10(pval),]
}

## Using bcftools

query_rsid_bcftools <- function(rsid, vcffile)
{
	bcftools <- get_bcftools_binary()
	tmp <- tempfile()
	write.table(unique(rsid), file=paste0(tmp, ".snplist"), row=F, col=F, qu=F)
	cmd <- sprintf("%s view -i'ID=@%s.snplist' %s > %s.vcf", bcftools, tmp, vcffile, tmp)
	system(cmd)
	o <- readVcf(paste0(tmp, ".vcf"))
	unlink(paste0(tmp, ".vcf"))
	unlink(paste0(tmp, ".snplist"))
	return(o)
}

query_pval_bcftools <- function(pval, vcffile, id=NULL)
{
	bcftools <- get_bcftools_binary()
	if(is.null(id))
	{
		id <- samples(scanVcfHeader(vcffile))
	}
	stopifnot(length(id) == 1)	
	tmp <- tempfile()
	write.table(unique(rsid), file=paste0(tmp, ".snplist"), row=F, col=F, qu=F)
	cmd <- sprintf("%s view -s %s -i 'FORMAT/LP > %s' %s > %s.vcf", bcftools, id, -log10(pval), vcffile, tmp)
	system(cmd)
	o <- readVcf(paste0(tmp, ".vcf"))
	unlink(paste0(tmp, ".vcf"))
	unlink(paste0(tmp, ".snplist"))
	return(o)
}

query_chrompos_bcftools <- function(chrompos, vcffile)
{
	bcftools <- get_bcftools_binary()
	chrompos <- parse_chrompos(chrompos)
	chrompos %>% as.data.frame
	tmp <- tempfile()
	write.table(as.data.frame(chrompos)[,1:3], file=paste0(tmp, ".snplist"), sep="\t", row=F, col=F, qu=F)

	cmd <- sprintf("%s view -R %s.snplist %s > %s.vcf", bcftools, tmp, vcffile, tmp)
	system(cmd)
	o <- readVcf(paste0(tmp, ".vcf"))
	unlink(paste0(tmp, ".vcf"))
	unlink(paste0(tmp, ".snplist"))
	return(o)
}


granges_from_df <- function(df)
{
	GRanges(seqnames=df$chrom, ranges=IRanges(start=df$start, end=df$end))
}


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


granges_from_vcf <- function(vcf, id=NULL)
{
	if(is.null(id))
	{
		id <- samples(header(vcf))
	}
	stopifnot(length(id) == 1)
	out <- VariantAnnotation::expand(vcf[,geno]) %>% 
	geno() %>%
	as.list() %>%
	lapply(., unlist) %>%
	bind_cols()
	a <- rowRanges(vcf)
	values(a) <- cbind(values(a), out)
	values(a)$id <- id

	return(a)
}


query_chrompos_file <- function(chrompos, vcffile, id=NULL, build="hg19")
{
	chrompos <- parse_chrompos(chrompos)
	tab <- TabixFile(vcffile)
	vcf <- readVcf(tab, build, param=chrompos)
	granges_from_vcf(vcf, id)
}


query_rsid_file <- function(rsid, vcffile)
{

}


query_pval_file <- function(pval, vcffile)
{

}


query_chrompos_vcf <- function(chrompos, vcffile)
{
	chrompos <- parse_chrompos(chrompos)
		
}


query_rsid_vcf <- function(rsid, vcffile)
{

}


query_pval_vcf <- function(pval, vcffile)
{

}







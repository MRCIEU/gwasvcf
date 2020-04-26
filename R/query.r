#' Query data from vcf file
#'
#' Read in GWAS summary data with filters on datasets (if multiple datasets per file) and/or chromosome/position, rsids or pvalues. Chooses most optimal choice for the detected operating system. Typically chrompos searches are the fastest. On Windows, rsid or pvalue filters from a file will be slow. 
#'
#' @param vcf VCF object e.g. output from VariantAnnotation::readVcf or gwasvcftools::query_vcf
#' @param chrompos Either vector of chromosome and position ranges e.g. "1:1000" or "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.
#' @param rsid Vector of rsids
#' @param pval  P-value threshold (NOT -log10)
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#' @param rsidx Path to rsidx index file
#' @param build ="GRCh37" Build of vcffile
#' @param os The operating system. Default is as detected. Determines the method used to perform query
#' @param proxies ="no" If SNPs are absent then look for proxies (yes) or not (no). Can also mask all target SNPs and only return proxies (only), for testing purposes. Currently only possible if querying rsid.
#' @param bfile =path to plink bed/bim/fam ld reference panel
#' @param dbfile =path to sqlite tag snp database
#' @param tag_kb =5000 Proxy parameter
#' @param tag_nsnp =5000 Proxy parameter
#' @param tag_r2 =0.6 Proxy parameter
#' @param threads =1 NUmber of threads
#'
#' @export
#' @return vcf object
query_gwas <- function(vcf, chrompos=NULL, rsid=NULL, pval=NULL, id=NULL, rsidx=NULL, build="GRCh37", os=Sys.info()[['sysname']], proxies="no", bfile=NULL, dbfile=NULL, tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1)
{
	if(is.character(vcf))
	{
		stopifnot(file.exists(vcf))
		fileflag <- TRUE
	} else {
		stopifnot(class(vcf) %in% c("CollapsedVCF", "ExpandedVCF"))
		fileflag <- FALSE
	}
	if(sum(c(!is.null(chrompos), !is.null(rsid), !is.null(pval))) != 1)
	{
		stop("Must specify filters only for one of chrompos, rsid or pval")
	}

	if(proxies != "no")
	{
		if(is.null(rsid))
		{
			stop("Proxies can only be searched for if rsid query specified")
		}
	}

	if(!is.null(chrompos))
	{
		if(!fileflag)
		{
			return(query_chrompos_vcf(chrompos, vcf))
		} else {
			if(!check_bcftools())
			{
				return(query_chrompos_file(chrompos, vcf, id, build))
			} else {
				return(query_chrompos_bcftools(chrompos, vcf, id))
			}
		}
	}

	if(!is.null(rsid))
	{
		stopifnot(proxies %in% c("yes", "no", "only"))
		if(proxies != "no")
		{
			return(proxy_match(vcf, rsid, bfile=bfile, dbfile=dbfile, proxies=proxies, tag_kb=tag_kb, tag_nsnp=tag_nsnp, tag_r2=tag_r2, threads=threads))
		}
		if(!fileflag)
		{
			return(query_rsid_vcf(rsid, vcf))
		} else {
			if(!is.null(rsidx))
			{
				return(query_rsid_rsidx(rsid, vcf, id, rsidx))
			}
			if(!check_bcftools())
			{
				return(query_rsid_file(rsid, vcf, id, build))
			} else {
				return(query_rsid_bcftools(rsid, vcf, id))
			}
		}
	}

	if(!is.null(pval))
	{
		if(!fileflag)
		{
			return(query_pval_vcf(pval, vcf, id))
		} else {
			if(!check_bcftools())
			{
				return(query_pval_file(pval, vcf, id, build))
			} else {
				return(query_pval_bcftools(pval, vcf, id))
			}
		}		
	}
}




df_to_granges <- function(df)
{
	GenomicRanges::GRanges(seqnames=df[["chrom"]], ranges=IRanges::IRanges(start=df[["start"]], end=df[["end"]]))
}






#' Parse chromosome:position
#'
#' Takes data frame or vector of chromosome position ranges and parses to granges object
#'
#' @param chrompos Either vector of chromosome and position ranges e.g. "1:1000" or "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.
#' @param radius Add radius to the specified positions. Default = NULL 
#'
#' @export
#' @return GRanges object
parse_chrompos <- function(chrompos, radius=NULL)
{

	if("GRanges" %in% class(chrompos))
	{
		if(!is.null(radius))
		{
			chrompos <- GenomicRanges::GRanges(
				seqnames = GenomeInfoDb::seqnames(chrompos),
				ranges = IRanges::IRanges(
					start = pmax(chrompos@start - radius, 0),
					end = chrompos@end + radius
				),
				strand = chrompos@strand
			)
		}
		return(chrompos)
	} else if(is.data.frame(chrompos)) {
		stopifnot(is.data.frame(chrompos))
		stopifnot(all(c("chrom", "start", "end") %in% names(chrompos)))
		return(df_to_granges(chrompos))
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
	if(!is.null(radius))
	{
		pos1 <- pmax(0, pos1 - radius)
		pos2 <- pos2 + radius
	}
	return(df_to_granges(data.frame(chrom, start=pos1, end=pos2, stringsAsFactors=FALSE)))
}




#' Query vcf file, extracting by chromosome and position
#'
#'
#' @param chrompos Either vector of chromosome and position ranges e.g. "1:1000" or "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.
#' @param vcffile Path to .vcf.gz GWAS summary data file
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#' @param build Default="GRCh37" Build of vcffile
#'
#' @export
#' @return VCF object
query_chrompos_file <- function(chrompos, vcffile, id=NULL, build="GRCh37")
{
	chrompos <- parse_chrompos(chrompos)
	if(!is.null(id))
	{
		param <- VariantAnnotation::ScanVcfParam(which=chrompos, samples=id)
	} else {
		param <- VariantAnnotation::ScanVcfParam(which=chrompos)
	}
	tab <- Rsamtools::TabixFile(vcffile)
	vcf <- VariantAnnotation::readVcf(tab, build, param=chrompos)
	return(vcf)
}


#' Query vcf file, extracting by rsid
#'
#' @param rsid Vector of rsids. Use DBSNP build (???)
#' @param vcffile Path to .vcf.gz GWAS summary data file
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#' @param build Default="GRCh37" Build of vcffile
#'
#' @export
#' @return VCF object
query_rsid_file <- function(rsid, vcffile, id=NULL, build="GRCh37")
{
	message("Note, this is much slower than searching by chromosome/position (e.g. see query_chrompos_file)")
	vcf <- Rsamtools::TabixFile(vcffile)
	fil <- function(x)
	{
		grepl(paste(rsid, collapse="|"), x)
	}

	tempfile <- tempfile()
	VariantAnnotation::filterVcf(vcf, build, tempfile, prefilters=S4Vectors::FilterRules(list(fil=fil)), verbose=TRUE)
	if(!is.null(id))
	{
		o <- VariantAnnotation::readVcf(tempfile, param=VariantAnnotation::ScanVcfParam(samples=id))
	} else {
		o <- VariantAnnotation::readVcf(tempfile)
	}
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
#' @param build Default="GRCh37"
#'
#' @export
#' @return VCF object
query_pval_file <- function(pval, vcffile, id=NULL, build="GRCh37")
{
	if(is.null(id))
	{
		id <- VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(vcffile))
	}
	stopifnot(length(id) == 1)
	message("Filtering p-value based on id ", id)
	message("Note, this is much slower than searching by chromosome/position (e.g. see query_chrompos_file)")
	vcf <- Rsamtools::TabixFile(vcffile)
	fil <- function(x)
	{
		VariantAnnotation::geno(x)[["LP"]][,id,drop=TRUE] > -log10(pval)
	}
	tempfile <- tempfile()
	VariantAnnotation::filterVcf(vcf, build, tempfile, filters=S4Vectors::FilterRules(list(fil=fil)), verbose=TRUE)
	o <- VariantAnnotation::readVcf(tempfile)
	unlink(tempfile)
	return(o)
}


#' Query chrompos from vcf object
#'
#' @param chrompos Either vector of chromosome and position ranges e.g. "1:1000" or "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.
#' @param vcf VCF object (e.g. from readVcf)
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#'
#' @export
#' @return VCF object
query_chrompos_vcf <- function(chrompos, vcf, id=NULL)
{
	if(is.null(id))
	{
		id <- VariantAnnotation::samples(VariantAnnotation::header(vcf))
	}
	colid <- which(VariantAnnotation::samples(VariantAnnotation::header(vcf)) == id)
	chrompos <- parse_chrompos(chrompos)
	i <- IRanges::findOverlaps(SummarizedExperiment::rowRanges(vcf), chrompos) %>% S4Vectors::queryHits() %>% unique %>% sort
	vcf[i,colid]
}


#' Query rsid from vcf object
#'
#' @param rsid Vector of rsids
#' @param vcf VCF object (e.g. from readVcf)
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#'
#' @export
#' @return VCF object
query_rsid_vcf <- function(rsid, vcf, id=NULL)
{
	if(is.null(id))
	{
		id <- VariantAnnotation::samples(VariantAnnotation::header(vcf))
	}
	colid <- which(VariantAnnotation::samples(VariantAnnotation::header(vcf)) == id)
	vcf[rownames(vcf) %in% rsid,colid]
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
		id <- VariantAnnotation::samples(VariantAnnotation::header(vcf))
	}
	stopifnot(length(id) == 1)
	colid <- which(VariantAnnotation::samples(VariantAnnotation::header(vcf)) == id)
	vcf[VariantAnnotation::geno(vcf)[["LP"]][,colid,drop=TRUE] > -log10(pval),colid]
}


#' Query 
#'
#' @param rsid Vector of rsids
#' @param vcffile Path to .vcf.gz GWAS summary data file
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#'
#' @export
#' @return VCF object
query_rsid_bcftools <- function(rsid, vcffile, id=NULL)
{
	stopifnot(check_bcftools())
	bcftools <- options()[["tools_bcftools"]]
	if(is.null(id))
	{
		id <- VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(vcffile))
	}
	id <- paste(id, collapse=",")
	tmp <- tempfile()
	utils::write.table(unique(rsid), file=paste0(tmp, ".snplist"), row=F, col=F, qu=F)
	cmd <- sprintf("%s view -s %s -i'ID=@%s.snplist' %s > %s.vcf", bcftools, id, tmp, vcffile, tmp)
	system(cmd)
	o <- VariantAnnotation::readVcf(paste0(tmp, ".vcf"))
	unlink(paste0(tmp, ".vcf"))
	unlink(paste0(tmp, ".snplist"))
	return(o)
}

#' Query p-value using bcftools
#'
#' @param pval P-value threshold (NOT -log10)
#' @param vcffile Path to .vcf.gz GWAS summary data file
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#'
#' @export
#' @return vcf object
query_pval_bcftools <- function(pval, vcffile, id=NULL)
{
	stopifnot(check_bcftools())
	bcftools <- options()[["tools_bcftools"]]
	if(is.null(id))
	{
		id <- VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(vcffile))
	}
	id <- paste(id, collapse=",")
	tmp <- tempfile()
	cmd <- sprintf("%s view -s %s -i 'FORMAT/LP > %s' %s > %s.vcf", bcftools, id, -log10(pval), vcffile, tmp)
	system(cmd)
	o <- VariantAnnotation::readVcf(paste0(tmp, ".vcf"))
	unlink(paste0(tmp, ".vcf"))
	return(o)
}

#' Query chromosome and position using bcftools
#'
#' @param chrompos Either vector of chromosome and position ranges e.g. "1:1000" or "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.
#' @param vcffile Path to .vcf.gz GWAS summary data file
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#'
#' @export
#' @return vcf object
query_chrompos_bcftools <- function(chrompos, vcffile, id=NULL)
{
	stopifnot(check_bcftools())
	bcftools <- options()[["tools_bcftools"]]
	if(is.null(id))
	{
		id <- VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(vcffile))
	}
	idclause <- ifelse(length(id) == 0, "", paste0("-s ", paste(id, collapse=",")))

	chrompos <- parse_chrompos(chrompos)
	chrompos %>% as.data.frame
	tmp <- tempfile()
	utils::write.table(as.data.frame(chrompos)[,1:3], file=paste0(tmp, ".snplist"), sep="\t", row=F, col=F, qu=F)

	cmd <- sprintf(paste0("%s view %s -R %s.snplist %s > %s.vcf"), bcftools, idclause, tmp, vcffile, tmp)
	system(cmd)
	o <- VariantAnnotation::readVcf(paste0(tmp, ".vcf"))
	unlink(paste0(tmp, ".vcf"))
	unlink(paste0(tmp, ".snplist"))
	return(o)
}


#' Query rsid from file using rsidx index
#'
#' See create_rsidx_index
#'
#' @param rsid Vector of rsids
#' @param vcffile Path to .vcf.gz GWAS summary data file
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#' @param rsidx Path to rsidx index file
#'
#' @export
#' @return vcf object
query_rsid_rsidx <- function(rsid, vcffile, id=NULL, rsidx)
{
	out <- query_rsidx(rsid, rsidx)
	return(
		query_gwas(vcffile, chrompos=data.frame(chrom=out$chrom, start=out$coord, end=out$coord), id=id)
	)
}

#' Query rsidx
#'
#' @param rsid Vector of rsids
#' @param rsidx Path to rsidx index file
#'
#' @export
#' @return data frame
query_rsidx <- function(rsid, rsidx)
{
	conn <- RSQLite::dbConnect(RSQLite::SQLite(), rsidx)
	numid <- gsub("rs", "", rsid) %>% paste(.data, collapse=",")
	query <- paste0("SELECT DISTINCT * FROM rsid_to_coord WHERE rsid IN (", numid, ")")
	out <- RSQLite::dbGetQuery(conn, query)
	RSQLite::dbDisconnect(conn)
	return(out)
}



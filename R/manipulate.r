#' Create GWAS vcf
#'
#' @param chrom chrom vector
#' @param pos pos vector
#' @param nea nea vector
#' @param ea ea vector
#' @param snp Optional vector
#' @param ea_af Optional vector
#' @param effect Optional vector
#' @param se Optional vector
#' @param pval Optional vector
#' @param n Optional vector
#' @param ncase Optional vector
#' @param name Optional vector
#'
#' @export
#' @return vcf object
create_vcf <- function(chrom, pos, nea, ea, snp=NULL, ea_af=NULL, effect=NULL, se=NULL, pval=NULL, n=NULL, ncase=NULL, name=NULL)
{
	stopifnot(length(chrom) == length(pos))
	if(is.null(snp))
	{
		snp <- paste0(chrom, ":", pos)
	}
	nsnp <- length(chrom)
	gen <- list()
	if(!is.null(ea_af)) gen[["AF"]] <- matrix(ea_af, nsnp)
	if(!is.null(effect)) gen[["ES"]] <- matrix(effect, nsnp)
	if(!is.null(se)) gen[["SE"]] <- matrix(se, nsnp)
	if(!is.null(pval)) gen[["LP"]] <- matrix(-log10(pval), nsnp)
	if(!is.null(n)) gen[["SS"]] <- matrix(n, nsnp)
	if(!is.null(ncase)) gen[["NC"]] <- matrix(ncase, nsnp)
	gen <- S4Vectors::SimpleList(gen)

	gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start=pos, end=pos + pmax(nchar(nea), nchar(ea)) - 1, names=snp))
	coldata <- S4Vectors::DataFrame(Samples = length(name), row.names=name)

	hdr <- VariantAnnotation::VCFHeader(
		header = IRanges::DataFrameList(
			fileformat = S4Vectors::DataFrame(Value="VCFv4.2", row.names="fileformat")
		),
		sample = name
	)
	VariantAnnotation::geno(hdr) <- S4Vectors::DataFrame(
		Number = c("A", "A", "A", "A", "A", "A"),
		Type = c("Float", "Float", "Float", "Float", "Float", "Float"),
		Description = c(
			"Effect size estimate relative to the alternative allele",
			"Standard error of effect size estimate",
			"-log10 p-value for effect estimate",
			"Alternate allele frequency in the association study",
			"Sample size used to estimate genetic effect",
			"Number of cases used to estimate genetic effect"
		),
		row.names=c("ES", "SE", "LP", "AF", "SS", "NC")
	)
	VariantAnnotation::geno(hdr) <- subset(VariantAnnotation::geno(hdr), rownames(VariantAnnotation::geno(hdr)) %in% names(gen))

	vcf <- VariantAnnotation::VCF(
		rowRanges = gr,
		colData = coldata,
		exptData = list(
			header = hdr
		),
		geno = gen
	)
	VariantAnnotation::alt(vcf) <- Biostrings::DNAStringSetList(as.list(ea))
	VariantAnnotation::ref(vcf) <- Biostrings::DNAStringSet(nea)
	VariantAnnotation::fixed(vcf)$FILTER <- "PASS"
	return(sort(vcf))
}

#' Merge two GWAS VCF objects
#'
#' Returns merged intersection of two VCF objects
#'
#' @param a VCF object
#' @param b VCF object
#'
#' @export
#' @return SimpleList of VCF objects
merge_vcf <- function(a, b)
{
	a <- VariantAnnotation::expand(a)
	b <- VariantAnnotation::expand(b)
	# o <- SummarizedExperiment::findOverlaps(a, b)
	o <- dplyr::tibble(
		from = which(names(a) %in% names(b)),
		to = match(names(a)[from], names(b))
	)
	a <- a[o[["from"]],]
	b <- b[o[["to"]],]
	allele_match <- VariantAnnotation::ref(a) == VariantAnnotation::ref(b) & VariantAnnotation::alt(a) == VariantAnnotation::alt(b)
	switch <- VariantAnnotation::ref(a) == VariantAnnotation::alt(b) & VariantAnnotation::ref(b) == VariantAnnotation::alt(a)
	if(any(switch))
	{
		for(i in 1:ncol(VariantAnnotation::geno(b)[["ES"]]))
		{
			VariantAnnotation::geno(b)[["ES"]][,i][switch] <- lapply(VariantAnnotation::geno(b)[["ES"]][,i][switch], function(x) x * -1)
		}
	}
	a <- a[allele_match | switch, ]
	b <- b[allele_match | switch, ]

	ab <- a
	temp <- lapply(names(VariantAnnotation::geno(ab)), function(x) rbind(VariantAnnotation::geno(ab)[x], VariantAnnotation::geno(b)[x])) %>% S4Vectors::SimpleList
	names(temp) <- names(VariantAnnotation::geno(ab))
	VariantAnnotation::geno(ab) <- temp

	h <- VariantAnnotation::header(a)
	VCFHeader(
		reference = VariantAnnotation::reference(h),
		samples = c(VariantAnnotation::samples(h), VariantAnnotation::samples(VariantAnnotation::header(b))),
		meta = VariantAnnotation::meta(h)
	)
	VariantAnnotation::samples(h) <- c(VariantAnnotation::samples(h), VariantAnnotation::samples(VariantAnnotation::header(b)))

	return(S4Vectors::SimpleList())
}



#' Convert vcf format to granges format
#'
#' @param vcf Output from readVcf
#' @param id Only accepts one ID, so specify here if there are multiple GWAS datasets in the vcf
#'
#' @export
#' @return GRanges object
vcf_to_granges <- function(vcf, id=NULL)
{
	stopifnot(class(vcf) %in% c("ExpandedVCF", "CollapsedVCF"))
	if(length(vcf) == 0)
	{
		message("VCF has length 0")
		return(NULL)
	}
	if(is.null(id))
	{
		id <- VariantAnnotation::samples(VariantAnnotation::header(vcf))
	}
	stopifnot(length(id) == 1)
	vcf <- VariantAnnotation::expand(vcf)
	a <- SummarizedExperiment::rowRanges(vcf)
	a$`REF` <- as.character(a$`REF`)
	a$`ALT` <- as.character(a$`ALT`)

	if(length(VariantAnnotation::geno(vcf)) == 0)
	{
		return(a)
	} else {
		out <- VariantAnnotation::expand(vcf) %>% 
			VariantAnnotation::geno() %>%
			as.list() %>%
			lapply(., function(x) unlist(x[,id,drop=TRUE])) %>%
			dplyr::bind_cols()
		S4Vectors::values(a) <- cbind(S4Vectors::values(a), out)
		S4Vectors::values(a)[["id"]] <- id

		if("TotalCases" %in% names(VariantAnnotation::meta(VariantAnnotation::header(vcf))$SAMPLE))
		{
			S4Vectors::values(a)[["NC"]] <- as.numeric(VariantAnnotation::meta(VariantAnnotation::header(vcf))$SAMPLE$TotalCases) %>% rep(., length(a))
			S4Vectors::values(a)[["SS"]] <- as.numeric(VariantAnnotation::meta(VariantAnnotation::header(vcf))$SAMPLE$TotalCases) + as.numeric(VariantAnnotation::meta(VariantAnnotation::header(vcf))$SAMPLE$TotalControls) %>% rep(., length(a))
		} else if("TotalControls" %in% names(VariantAnnotation::meta(VariantAnnotation::header(vcf))$SAMPLE)) {
			S4Vectors::values(a)[["SS"]] <- as.numeric(VariantAnnotation::meta(VariantAnnotation::header(vcf))$SAMPLE$TotalControls) %>% rep(., length(a))
		}
		return(a)
	}
}


#' Convert vcf format to tibble (data frame)
#'
#' @param vcf Output from readVcf
#' @param id Only accepts one ID, so specify here if there are multiple GWAS datasets in the vcf
#'
#' @export
#' @return GRanges object
vcf_to_tibble <- function(vcf, id=NULL)
{
	a <- vcf_to_granges(vcf, id)
	if(is.null(a))
	{
		return(dplyr::tibble())
	}
	S4Vectors::values(a)[["rsid"]] <- names(a)
	return(dplyr::as_tibble(a))
}


#' Reduce list of VCFs to intersecting regions
#'
#' @param vcflist List of VCF objects, or list of VCF filenames, or mix of VCF objects and filenames
#' @param chrompos Either vector of chromosome and position ranges e.g. "1:1000" or "1:1000-2000", or data frame with columns `chrom`, `start`, `end`. 
#'
#' @export
#' @return List of VCFs
vcflist_overlaps <- function(vcflist, chrompos)
{
	stopifnot(is.list(vcflist))
	if(!is.null(chrompos))
	{
		chrompos <- parse_chrompos(chrompos)
		vcflist <- lapply(vcflist, function(x)
		{
			query_gwas(x, chrompos)
		})
	} else {
		vcflist <- lapply(1:length(vcflist), function(i)
		{
			x <- vcflist[[i]]
			if(class(x) %in% c("CollapsedVCF", "ExpandedVCF"))
			{
				if(is.null(chrompos))
				{
					return(x)
				} else {
					return(query_gwas(x, chrompos))
				}
			}
			if(class(x) == "character")
			{
				if(is.null(chrompos))
				{
					return(VariantAnnotation::readVcf(x))
				} else {
					return(query_gwas(x, chrompos))
				}
			}
			stop("Item ", i, " in vcflist is neither VCF object nor path to VCF file")
		})
	}

	# collapse indels for sorting purposes
	vcflist <- lapply(vcflist, function(x)
	{
		SummarizedExperiment::end(x) <- SummarizedExperiment::start(x)
	
		# Simple approach to avoid duplicate positions due to snps and indels
		x <- x[!duplicated(SummarizedExperiment::start(x))]
		return(x)
	})



	o <- Reduce(IRanges::subsetByOverlaps, lapply(vcflist, SummarizedExperiment::rowRanges))
	vcflist <- lapply(vcflist, function(x) IRanges::subsetByOverlaps(x, o))
	return(vcflist)
}

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
	return(vcf)
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
	a <- expand(a)
	b <- expand(b)
	# o <- SummarizedExperiment::findOverlaps(a, b)
	o <- dplyr::tibble(
		from = which(names(a) %in% names(b)),
		to = match(names(a)[from], names(b))
	)
	a <- a[o[["from"]],]
	b <- b[o[["to"]],]
	allele_match <- ref(a) == ref(b) & alt(a) == alt(b)
	switch <- ref(a) == alt(b) & ref(b) == alt(a)
	if(any(switch))
	{
		for(i in 1:ncol(geno(b)[["ES"]]))
		{
			geno(b)[["ES"]][,i][switch] <- lapply(geno(b)[["ES"]][,i][switch], function(x) x * -1)
		}
	}
	a <- a[allele_match | switch, ]
	b <- b[allele_match | switch, ]

	ab <- a
	temp <- lapply(names(geno(ab)), function(x) rbind(geno(ab)[x], geno(b)[x])) %>% SimpleList
	names(temp) <- names(geno(ab))
	geno(ab) <- temp

	h <- header(a)
	VCFHeader(
		reference = reference(h),
		samples = c(samples(h), samples(header(b))),
		meta = meta(h)
	)
	samples(h) <- c(samples(h), samples(header(b)))

	return(S4Vectors::SimpleList())
}



#' Create exposure or outcome data format for TwoSampleMR from vcf
#'
#' @param vcf VCF object
#' @param type ="exposure" or "outcome"
#'
#' @export
#' @return data frame
vcf_to_TwoSampleMR <- function(vcf, type="exposure")
{
	a <- vcf %>% vcf_to_granges
	a[["SNP"]] <- names(a)
	a <- dplyr::as_tibble(a)
	if(!"ES" %in% names(a)) a[["ES"]] <- NA
	if(!"SE" %in% names(a)) a[["SE"]] <- NA
	if(!"LP" %in% names(a)) a[["LP"]] <- NA
	if(!"SS" %in% names(a)) a[["SS"]] <- NA
	if(!"NC" %in% names(a)) a[["NC"]] <- NA
	if(!"id" %in% names(a)) a[["id"]] <- NA
	a[["LP"]] <- 10^-a[["LP"]]
	a[["NCONT"]] <- a[["SS"]] - a[["NC"]]
	TwoSampleMR::format_data(
		a, type=type,
		snp_col="SNP",
		effect_allele_col="ALT",
		other_allele_col="REF",
		eaf_col="AF",
		chr_col="CHROM",
		pos_col="POS",
		beta_col="ES",
		se_col="SE",
		pval_col="LP",
		samplesize_col="SS",
		ncase_col="NC",
		ncontrol_col="NCONT",
		phenotype_col="id"
	)
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
	if(is.null(id))
	{
		id <- VariantAnnotation::samples(VariantAnnotation::header(vcf))
	}
	stopifnot(length(id) == 1)
	a <- SummarizedExperiment::rowRanges(vcf)
	a$`ALT` <- unlist(a$`ALT`)

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

		if("TotalCases" %in% names(meta(header(vcf))$SAMPLE))
		{
			a[["NC"]] <- as.numeric(meta(header(vcf))$SAMPLE$TotalCases) %>% rep(., nrow(a))
			a[["SS"]] <- as.numeric(meta(header(vcf))$SAMPLE$TotalCases) + as.numeric(meta(header(vcf))$SAMPLE$TotalControls) %>% rep(., nrow(a))
		} else if("TotalControls" %in% names(meta(header(vcf))$SAMPLE)) {
			a[["SS"]] <- as.numeric(meta(header(vcf))$SAMPLE$TotalControls) %>% rep(., nrow(a))
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
	a[["SNP"]] <- names(vcf)
	return(dplyr::as_tibble(a))
}


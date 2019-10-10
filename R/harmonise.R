#' Read in GWAS dataset
#'
#'
#' @param filename <what param does>
#' @param skip <what param does>
#' @param delimiter <what param does>
#' @param gzipped <what param does>
#' @param snp <what param does>
#' @param nea <what param does>
#' @param ea <what param does>
#' @param ea_af <what param does>
#' @param effect <what param does>
#' @param se <what param does>
#' @param pval <what param does>
#' @param n <what param does>
#' @param info <what param does>
#' @param z <what param does>
#'
#' @export
#' @return data frame with log attributes
read_gwas <- function(filename, skip, delimiter, gzipped, snp, nea, ea, ea_af, effect, se, pval, n, info, z)
{
	if(gzipped)
	{
		# dat <- data.table::fread(paste0("gunzip -c ", filename), header=FALSE, skip=skip, sep=delimiter)
		dat <- data.table::fread(filename, header=FALSE, skip=skip, sep=delimiter)
	} else {
		dat <- data.table::fread(filename, header=FALSE, skip=skip, sep=delimiter)
	}
	nc <- ncol(dat)
	jlog <- list()
	jlog[['total_variants']] <- nrow(dat)

	if(snp == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		snp <- ncol(dat)
	}
	if(nea == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		nea <- ncol(dat)
	}
	if(ea == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		ea <- ncol(dat)
	}
	if(ea_af == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		ea_af <- ncol(dat)
	}
	if(effect == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		effect <- ncol(dat)
	}
	if(se == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		se <- ncol(dat)
	}
	if(pval == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		pval <- ncol(dat)
	}
	if(n == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		n <- ncol(dat)
	}
	if(info == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		info <- ncol(dat)
	}
	if(z == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		z <- ncol(dat)
	}
	o <- TwoSampleMR::format_data(
		dat, 
		type="outcome", 
		phenotype_col="outcome",
		snp_col=names(dat)[snp],
		beta_col=names(dat)[effect],
		se_col=names(dat)[se],
		effect_allele_col=names(dat)[ea],
		other_allele_col=names(dat)[nea],
		eaf_col=names(dat)[ea_af],
		pval_col=names(dat)[pval],
		samplesize_col=names(dat)[n],
		info_col=names(dat)[info],
		z_col=names(dat)[z],
	)
	print(head(o))
	o$beta.outcome[!is.finite(o$beta.outcome)] <- NA
	o$se.outcome[!is.finite(o$se.outcome)] <- NA
	o$pval.outcome[!is.finite(o$pval.outcome)] <- NA
	o$eaf.outcome[!is.finite(o$eaf.outcome)] <- NA
	o$info.outcome[!is.finite(o$info.outcome)] <- NA
	o$z.outcome[!is.finite(o$z.outcome)] <- NA
	if(all(is.na(o$info.outcome)))
	{
		o <- subset(o, select=-c(info.outcome))
	}
	if(all(is.na(o$z.outcome)))
	{
		o <- subset(o, select=-c(z.outcome))
	}
	jlog[['variants_not_read']] <- nrow(dat) - nrow(o)
	ind <- is.finite(o$beta.outcome) & is.finite(o$pval.outcome)
	o <- o[ind,]
	jlog[['variants_with_missing_stats']] <- sum(!ind)
	jlog[['variants_with_missing_pvals']] <- sum(is.na(o$pval.outcome))
	o <- subset(o, !is.na(pval.outcome))
	attr(o, "log") <- jlog
	return(o)
}



#' Read in reference dataset
#'
#' @param reference_file <what param does>
#' @param snplist <what param does>
#' @param outfile <what param does>
#' @param remove_dup_rsids=TRUE <what param does>
#'
#' @export
#' @return data frame
read_reference <- function(reference_file, rsid=NULL, chrompos=NULL, remove_dup_rsids=TRUE)
{
	if(!is.null(rsid))
	{
		a <- query_gwas(reference_file, rsid=rsid)
	} else if(!is.null(chrompos)) {
		a <- query_gwas(reference_file, chrompos=chrompos)
	} else {
		a <- query_gwas(reference_file, chrompos=chrompos)
	}
	if(remove_dup_rsids)
	{
		a <- a[!duplicated(names(a)), ]
	}
	return(format_from_vcf(a))
}

#' Harmonise gwas alleles to be same as reference
#'
#' @param gwas <what param does>
#' @param reference <what param does>
#'
#' @export
#' @return data frame with attributes
harmonise_against_ref <- function(gwas, reference)
{
	# Check strand
	action <- is_forward_strand(gwas$SNP, gwas$effect_allele.outcome, gwas$other_allele.outcome, reference$SNP, reference$other_allele.exposure, reference$effect_allele.exposure, threshold=0.9)

	# Harmonise
	dat <- TwoSampleMR::harmonise_data(reference, gwas, action)
	jlog <- c(
		attr(gwas, "log"),
		as.list(attr(dat, "log"))
	)
	jlog[['harmonised_variants']] <- sum(dat$mr_keep)
	jlog[['variants_not_harmonised']] <- sum(!dat$mr_keep)

	cols <- c("SNP"="ID", "effect_allele.exposure"="ALT", "other_allele.exposure"="REF", "beta.outcome"="BETA", "se.outcome"="SE", "pval.outcome"="PVALUE", "eaf.outcome"="AF", "samplesize.outcome"="N", "z.outcome"="ZVALUE", "info.outcome"="INFO")
	if(! "z.outcome" %in% names(dat)) dat$z.outcome <- NA
	if(! "info.outcome" %in% names(dat)) dat$info.outcome <- NA

	dat <- subset(dat, select=names(cols))
	names(dat) <- cols

	names(ref)[names(ref) == "chr.exposure"] <- "CHROM"
	names(ref)[names(ref) == "pos.exposure"] <- "POS"

	dat <- dat %>%
		dplyr::inner_join(subset(ref, select=c(SNP,other_allele.exposure,effect_allele.exposure,CHROM,POS)), by=c("ID"="SNP", "REF"="other_allele.exposure", "ALT"="effect_allele.exposure"))
	jlog[['total_remaining_variants']] <- nrow(dat)
	attr(dat, "log") <- jlog
	return(dat)
}


#' Create exposure or outcome data format for TwoSampleMR from vcf
#'
#' @param vcf VCF object
#' @param type="exposure" or "outcome"
#'
#' @export
#' @return
format_from_vcf <- function(vcf, type="exposure")
{
	a <- vcf %>% vcf_to_granges
	a$SNP <- names(a)
	a <- as_tibble(a)
	if(!"ES" %in% names(a)) a$ES <- NA
	if(!"SE" %in% names(a)) a$SE <- NA
	if(!"LP" %in% names(a)) a$LP <- NA
	if(!"SS" %in% names(a)) a$SS <- NA
	if(!"NC" %in% names(a)) a$NC <- NA
	if(!"id" %in% names(a)) a$id <- NA
	a$LP <- 10^-a$LP
	a$NCONT <- a$SS - a$NC
	format_data(
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



#' Create format for HPC pipeline 
#'
#' Takes raw files and aligns them to reference. Important if files don't have chr:pos already
#'
#' @param harmonised Output from /code{harmonise_against_ref}
#' @param path Path to write out json file and txt file
#'
#' @export
#' @return NULL
write_out <- function(harmonised, path)
{
	j <- list(
		chr_col = 10,
		pos_col = 11,
		snp_col = 0,
		ea_col = 2,
		oa_col = 1,
		beta_col = 3,
		se_col = 4,
		ncontrol_col = 7,
		pval_col = 5,
		eaf_col = 6,
		delimiter = " ",
		header = TRUE,
		build = "GRCh37"
	)
	if(!all(is.na(harmonised$ZVALUE))) j$imp_z_col <- 8
	if(!all(is.na(harmonised$INFO))) j$imp_z_col <- 9
	if(!all(is.na(harmonised$NCASE))) j$ncase_col <- which(names(harmonised) == "NCASE")

	write_json(j, paste0(path, ".json"), auto_unbox=TRUE)
	gz1 <- gzfile(paste0(path, ".txt.gz"), "w")
	write.table(harmonised, gz1, row=FALSE, col=TRUE, qu=FALSE)
	close(gz1)

}


#' Check a GWAS dataset against a reference known to be on the forward strand
#'
#' 
#' Assuming reference data is all on forward strand, check if 
#' the GWAS is also.
#' Use some threshold e.g. if more than 90% of alleles don't 
#' need to be flipped then it's likely that the dataset is on
#' the forward strand
#'
#' This function can be used to evaluate how strict harmonisation should be
#' The trade off if you assume we are not on the forward strand then palindromic SNPs are dropped within a particular frequency range
#' But you could instead have some small probability of error for whether palindromic SNPs are on the forward strand, and avoid dropping too many variants.
#'
#' @param gwas_snp Vector of SNP names for the dataset being checked
#' @param gwas_a1 Vector of alleles
#' @param gwas_a2 Vector of alleles
#' @param ref_snp Vector of SNP names for the reference dataset
#' @param ref_a1 Vector of alleles
#' @param ref_a2 Vector of alleles
#' @param threshold=0.9 If the proportion of allele strands match is above this threshold, then declare the dataset to be on the forward strand
#'
#' @export
#' @return 1 = Forward strand; 2 = Not on forward strand
is_forward_strand <- function(gwas_snp, gwas_a1, gwas_a2, ref_snp, ref_a1, ref_a2, threshold=0.9)
{
	requireNamespace("dplyr", quietly=TRUE)
	if(is.null(gwas_a1) | is.null(gwas_a2))
	{
		message("No info for both GWAS alleles")
		return(2)
	}

	if(1-(sum(is.na(gwas_a1)) / length(gwas_a1)) < threshold)
	{
		message("Too many missing values for gwas A1")
		return(2)
	}
	if(1-(sum(is.na(gwas_a2)) / length(gwas_a2)) < threshold)
	{
		message("Too many missing values for gwas A2")
		return(2)
	}

	g <- dplyr::data_frame(SNP=gwas_snp, A1=toupper(gwas_a1), A2=toupper(gwas_a2)) %>%
    		subset(!is.na(SNP) & !is.na(A1) & !is.na(A2))
	r <- dplyr::data_frame(SNP=ref_snp, A1=toupper(ref_a1), A2=toupper(ref_a2))

	gr <- dplyr::inner_join(g,r,by="SNP")
	index <- (gr$A1.x == gr$A1.y & gr$A2.x == gr$A2.y) | (gr$A1.x == gr$A2.y & gr$A2.x == gr$A1.y)
	diff <- gr[!index,]
	diff$A1.x[diff$A1.x == "C"] <- "g"
	diff$A1.x[diff$A1.x == "G"] <- "c"
	diff$A1.x[diff$A1.x == "T"] <- "a"
	diff$A1.x[diff$A1.x == "A"] <- "t"
	diff$A2.x[diff$A2.x == "C"] <- "g"
	diff$A2.x[diff$A2.x == "G"] <- "c"
	diff$A2.x[diff$A2.x == "T"] <- "a"
	diff$A2.x[diff$A2.x == "A"] <- "t"
	diff$A1.x <- toupper(diff$A1.x)
	diff$A2.x <- toupper(diff$A2.x)

	index2 <- (diff$A1.x == diff$A1.y & diff$A2.x == diff$A2.y) | (diff$A1.x == diff$A2.y & diff$A2.x == diff$A1.y)

	# Number that match initially
	message("SNPs that match: ", sum(index))
	message("SNPs that match after flipping: ", sum(index2))
	message("SNPs that never match: ", sum(!index2))

	prop <- 1 - sum(index2) / sum(index)
	message("Proportion on forward strand: ", prop)

	return(ifelse(prop > threshold, 1, 2))
}



check_null <- function(x, n)
{
	if(is.null(x))
	{
		return(rep(NA))
	}
}

#' <brief desc>
#'
#' <full description>
#'
#' @param dat <what param does>
#' @param snp_col <what param does>
#' @param chrom_col <what param does>
#' @param pos_col <what param does>
#' @param nea_col <what param does>
#' @param ea_col <what param does>
#' @param ea_af_col <what param does>
#' @param effect_col <what param does>
#' @param se_col <what param does>
#' @param pval_col <what param does>
#' @param n_col <what param does>
#' @param ncase_col <what param does>
#' @param info_col <what param does>
#' @param z_col <what param does>
#'
#' @export
#' @return
create_vcf <- function(chrom, pos, nea, ea, snp=NULL, ea_af=NULL, effect=NULL, se=NULL, pval=NULL, n=NULL, ncase=NULL, name=NULL)
{
	stopifnot(length(chrom) == length(pos))
	if(is.null(snp))
	{
		snp <- paste0(chrom, ":", pos)
	}
	nsnp <- length(chrom)
	gen <- list()
	if(!is.null(ea_af)) gen$AF <- matrix(ea_af, nsnp)
	if(!is.null(effect)) gen$ES <- matrix(effect, nsnp)
	if(!is.null(se)) gen$SE <- matrix(se, nsnp)
	if(!is.null(pval)) gen$LP <- matrix(-log10(pval), nsnp)
	if(!is.null(n)) gen$SS <- matrix(n, nsnp)
	if(!is.null(ncase)) gen$NC <- matrix(ncase, nsnp)
	gen <- SimpleList(gen)

	gr <- GRanges(chrom, IRanges(start=pos, end=pos + pmax(nchar(nea), nchar(ea)) - 1, names=snp))
	coldata <- S4Vectors::DataFrame(Samples = length(name), row.names=name)

	header <- VCFHeader(
		header = DataFrameList(
			fileformat = DataFrame(Value="VCFv4.2", row.names="fileformat")
		),
		sample = name
	)
	geno(header) <- DataFrame(
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
	) %>% subset(., rownames(.) %in% names(gen))

	vcf <- VCF(
		rowRanges = gr,
		colData = coldata,
		exptData = list(
			header = header,
			fixed = DataFrame(REF=DNAStringSet(nea), ALT=DNAStringSetList(as.list(ea)), QUAL=as.numeric(NA), FILTER="PASS")
		),
		geno = gen
	)
	return(vcf)
}


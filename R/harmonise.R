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
read_reference <- function(reference_file, snplist, outfile, remove_dup_rsids=TRUE)
{
	stopifnot(check_bcftools() == 0)
	out1 <- paste0(outfile, ".snplist")
	out2 <- paste0(outfile, ".ref")
	write.table(snplist, file=out1, row=F, col=F, qu=F)

	# Extract relevent info
	cmd <- paste0("time bcftools view -i'ID=@", out1, "' ", reference_file, " | bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' > ", out2)
	print(cmd)
	system(cmd)
	unlink(out1)
	a <- data.table::fread(out2, header=FALSE, sep="\t")
	unlink(out2)
	names(a) <- c("CHROM", "POS", "ID", "REF", "ALT", "AF")

	if(remove_dup_rsids)
	{
		a <- subset(a, !duplicated(ID))
	}
	a$beta <- 1
	a$se <- 0.1
	a$pval <- 0.1
	a <- TwoSampleMR::format_data(
		a,
		type="exposure",
		snp_col="ID",
		effect_allele_col="ALT",
		other_allele_col="REF",
		eaf_col="AF",
		chr_col="CHROM",
		pos_col="POS"
	)
	return(a)
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
	action <- gwasvcftools::is_forward_strand(gwas$SNP, gwas$effect_allele.outcome, gwas$other_allele.outcome, reference$SNP, reference$other_allele.exposure, reference$effect_allele.exposure, threshold=0.9)

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



#' Create format for HPC pipeline 
#'
#' Takes raw files and aligns them to reference. Important if files don't have chr:pos already
#'
#' @param harmonsied Output from /code{harmonise_against_ref}
#' @param path Path to write out json file and txt file
#'
#' @export
#' @return NULL
write_out <- function(harmonsied, path)
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

	write_json(j, paste0(path, ".json"))
	write.table(harmonised, file=paste0(path, ".txt"), row=FALSE, col=TRUE, qu=FALSE)
}

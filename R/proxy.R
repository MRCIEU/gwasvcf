#' Find LD proxies for a set of SNPs
#'
#' @param rsid list of rs IDs
#' @param bfile ld reference panel
#' @param tag_kb =5000 Proxy parameter
#' @param tag_nsnp =5000 Proxy parameter
#' @param tag_r2 =0.6 Proxy parameter
#' @param searchspace Optional list of rs IDs to use as potential proxies 
#' @param threads Number of threads to use (=1)
#' @param out temporary output file
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @return data frame
get_ld_proxies <- function(rsid, bfile, searchspace=NULL, tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1, out=tempfile())
{
	stopifnot(check_plink())
	searchspacename <- paste0(out, ".searchspace")
	targetsname <- paste0(out, ".targets")
	outname <- paste0(out, ".targets.ld.gz")
	utils::write.table(rsid, file=targetsname, row=FALSE, col=FALSE, qu=FALSE)
	if(!is.null(searchspace))
	{
		stopifnot(is.character(searchspace))
		utils::write.table(c(rsid, searchspace), file=searchspacename, row=F, col=F, qu=F)
		extract_param <- paste0(" --extract ", searchspacename)
	} else {
		extract_param <- " "
	}
	cmd <- paste0(options()[["tools_plink"]],
		" --bfile ", bfile, 
		extract_param,
		" --keep-allele-order ",
		" --r in-phase with-freqs gz",
		" --ld-snp-list ", targetsname,
		" --ld-window-kb ", tag_kb,
		" --ld-window ", tag_nsnp,
		" --out ", targetsname,
		" --threads ", threads,
		" 2>&1 > /dev/null"
	)
	system(cmd)

	ld <- data.table::fread(paste0("gunzip -c ", outname), header=TRUE) %>%
		dplyr::filter(`R`^2 > `tag_r2`) %>%
		dplyr::filter(`SNP_A` != `SNP_B`) %>%
		dplyr::mutate(PHASE=gsub("/", "", `PHASE`)) %>%
		subset(., nchar(`PHASE`) == 4)
	if(nrow(ld) == 0)
	{
		return(ld)
	}
	temp <- do.call(rbind, strsplit(ld[["PHASE"]], "")) %>% dplyr::as_tibble(., .name_repair="minimal")
	names(temp) <- c("A1", "B1", "A2", "B2")
	ld <- cbind(ld, temp) %>% dplyr::as_tibble(., .name_repair="minimal")
	# ld <- dplyr::arrange(ld, desc(abs(R))) %>%
	# 	dplyr::filter(!duplicated(SNP_A))
	ld <- dplyr::arrange(ld, dplyr::desc(abs(`R`)))
	unlink(searchspacename)
	unlink(targetsname)
	unlink(paste0(targetsname, c(".log", ".nosex")))
	unlink(outname)
	return(ld)
}


#' Extract SNPs from vcf file
#'
#' Finds proxies if necessary
#'
#' @param vcf vcf file name
#' @param rsid list of rs IDs
#' @param bfile ld reference panel
#' @param proxies ="yes" If SNPs are absent then look for proxies (yes) or not (no). Can also mask all target SNPs and only return proxies (only), for testing purposes
#' @param tag_kb =5000 Proxy parameter
#' @param tag_nsnp =5000 Proxy parameter
#' @param tag_r2 =0.6 Proxy parameter
#' @param threads Number of threads to use (=1)

#'
#' @export
#' @return data frame
proxy_match <- function(vcf, rsid, bfile, proxies="yes", tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1)
{

	os <- Sys.info()[['sysname']]
	if(proxies=="yes")
	{
		o <- query_gwas(vcf, rsid=rsid)
		missing <- rsid[!rsid %in% names(o)]
		if(length(missing) != 0)
		{
			searchspacename <- tempfile()
			if(is.character(vcf))
			{
				if(check_bcftools())
				{
					cmd <- paste0(options()[["tools_bcftools"]], " query -f'%ID\n' ", vcf, " > ", searchspacename)
					system(cmd)
					searchspace <- scan(searchspacename, what="character")
				} else {
					searchspace <- NULL
				}				
			} else {
				searchspace <- names(SummarizedExperiment::rowRanges(vcf))
			}
			ld <- get_ld_proxies(missing, bfile, searchspace=searchspace, tag_kb=tag_kb, tag_nsnp=tag_nsnp, tag_r2=tag_r2, threads=threads)
			if(nrow(ld) == 0)
			{
				return(o)
			}
		} else {
			return(o)
		}
	} else if(proxies == "only") {
			searchspacename <- tempfile()
			if(is.character(vcf))
			{
				if(check_bcftools())
				{
					cmd <- paste0(options()[["tools_bcftools"]], " query -f'%ID\n' ", vcf, " > ", searchspacename)
					system(cmd)
					searchspace <- scan(searchspacename, what="character")
				} else {
					searchspace <- NULL
				}				
			} else {
				searchspace <- names(SummarizedExperiment::rowRanges(vcf))
			}
		ld <- get_ld_proxies(rsid, bfile, searchspace=searchspace, tag_kb=tag_kb, tag_nsnp=tag_nsnp, tag_r2=tag_r2, threads=threads)
		if(nrow(ld) == 0)
		{
			return(VCF())
		}
	} else if(proxies == "no") {
		o <- query_gwas(vcf, rsid=rsid)
		return(o)
	} else {
		stop('proxies must be "yes", "no" or "only"')
	}
	if(!is.null(searchspace))
	{
		ld <- ld %>% dplyr::filter(!duplicated(SNP_A))
	}
	e <- query_gwas(vcf, rsid=ld[["SNP_B"]])

	if(is.null(searchspace))
	{
		ld <- subset(ld, `SNP_B` %in% names(e)) %>%
			dplyr::filter(!duplicated(`SNP_A`))		
	}
	e <- e[names(e) %in% ld[["SNP_B"]], ]
	index <- match(names(e), ld[["SNP_B"]])
	ld <- ld[index,]
	stopifnot(all(ld[["SNP_B"]] == names(e)))
	sign_index <- SummarizedExperiment::rowRanges(e)$`REF` == ld[["B1"]]

	gr <- GenomicRanges::GRanges(ld[["CHR_A"]], IRanges::IRanges(start=ld[["BP_A"]], end=ld[["BP_A"]], names=ld[["SNP_A"]]))
	fixeddat <- S4Vectors::DataFrame(
		paramRangeID=as.factor(rep(NA, nrow(ld))),  
		REF=Biostrings::DNAStringSet(ld[["A1"]]), 
		ALT=Biostrings::DNAStringSetList(as.list(ld[["A2"]])), 
		QUAL=as.numeric(NA), 
		FILTER="PASS"
	)
	prox <- VariantAnnotation::VCF(
		rowRanges = gr,
		colData = SummarizedExperiment::colData(e),
		info = VariantAnnotation::info(e),
		exptData = list(
			header = VariantAnnotation::header(e), 
			fixed = fixeddat
		),
		geno = S4Vectors::SimpleList(
			lapply(VariantAnnotation::geno(e), `dimnames<-`, NULL)
		)
	)

	VariantAnnotation::geno(VariantAnnotation::header(prox)) <- rbind(VariantAnnotation::geno(VariantAnnotation::header(prox)), 
		S4Vectors::DataFrame(row.names="PR", Number="1", Type="String", Description="Proxy rsid")
	)
	VariantAnnotation::geno(prox)[["ES"]][!sign_index] <- {unlist(VariantAnnotation::geno(prox)[["ES"]][!sign_index]) * -1} %>% as.list
	VariantAnnotation::geno(prox)[["PR"]] <- matrix(ld[["SNP_B"]], length(ld[["SNP_B"]]), 1)

	if(proxies == "only")
	{
		return(prox)
	} else {
		VariantAnnotation::geno(VariantAnnotation::header(o)) <- rbind(VariantAnnotation::geno(VariantAnnotation::header(o)), S4Vectors::DataFrame(row.names="PR", Number="1", Type="String", Description="Proxy rsid"))
		VariantAnnotation::geno(o)[["PR"]] <- matrix(rep(NA, length(o)), length(o), 1)
		return(BiocGenerics::rbind(o, prox))
	}
}


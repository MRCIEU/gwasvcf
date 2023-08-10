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
#' @importFrom rlang .data
#'
#' @export
#' @return data frame
get_ld_proxies <- function(rsid, bfile, searchspace=NULL, tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1, out=tempfile())
{
	stopifnot(check_plink())
	searchspacename <- paste0(out, ".searchspace")
	targetsname <- paste0(out, ".targets")
	outname <- paste0(out, ".targets.ld.gz")
	utils::write.table(rsid, file=targetsname, row.names = FALSE, col.names = FALSE, quote = FALSE)
	if(!is.null(searchspace))
	{
		stopifnot(is.character(searchspace))

		utils::write.table(unique(c(rsid, searchspace)), file=searchspacename, row.names = FALSE, col.names = FALSE, quote = FALSE)
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
		" --ld-window-r2 ", tag_r2,
		" --ld-window ", tag_nsnp,
		" --out ", targetsname,
		" --threads ", threads,
		" 2>&1 > /dev/null"
	)
	message("Finding proxies...")
	system(cmd)

	if (Sys.info()["sysname"] == "Windows") {
	  stop("Currently, this function only works on macOS and Linux")
	}
	if (!file.exists(outname)) {
		ld <- data.frame(CHR_A = integer(), BP_A = integer(), SNP_A = character(), MAF_A = double(), CHR_B = integer(), BP_B = integer(), 
		SNP_B = character(), PHASE = character(), MAF_B = double(), R = double())
		message("Index SNP not found in the reference panel")
		return(ld)
	}
	ld <- data.table::fread(paste0("gunzip -c ", outname), header=TRUE) %>%
		dplyr::as_tibble(.name_repair="minimal") %>%
		dplyr::filter(.data[["R"]]^2 > tag_r2) %>%
		dplyr::filter(.data[["SNP_A"]] != .data[["SNP_B"]]) %>%
		dplyr::mutate(PHASE=gsub("/", "", .data[["PHASE"]])) %>%
		dplyr::filter(nchar(.data[["PHASE"]]) == 4)
	unlink(searchspacename)
	unlink(targetsname)
	unlink(paste0(targetsname, c(".log", ".nosex")))
	unlink(outname)
	if(nrow(ld) == 0)
	{
		message("No proxies found")
		return(ld)
	}
	temp <- do.call(rbind, strsplit(ld[["PHASE"]], "")) %>% dplyr::as_tibble(.name_repair="minimal")
	names(temp) <- c("A1", "B1", "A2", "B2")
	ld <- cbind(ld, temp) %>% dplyr::as_tibble(.name_repair="minimal")
	# ld <- dplyr::arrange(ld, desc(abs(R))) %>%
	# 	dplyr::filter(!duplicated(SNP_A))
	ld <- dplyr::arrange(ld, dplyr::desc(abs(.data[["R"]])))
	message("Found ", nrow(ld), " proxies")
	return(ld)
}



#' Lookup LD proxies from sqlite database
#'
#' @param rsids List of rsids
#' @param dbfile path to dbfile
#' @param tag_r2 minimum r2 value
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
#' @return data frame
sqlite_ld_proxies <- function(rsids, dbfile, tag_r2)
{
	conn <- RSQLite::dbConnect(RSQLite::SQLite(), dbfile)
	numid <- gsub("rs", "", rsids) %>% paste(collapse=",")
	query <- paste0("SELECT DISTINCT * FROM tags WHERE SNP_A IN (", numid, ")")
	ld <- RSQLite::dbGetQuery(conn, query) %>% 
		dplyr::as_tibble(.name_repair="minimal") %>%
		dplyr::filter(.data[["R"]]^2 > tag_r2) %>%
		dplyr::filter(.data[["SNP_A"]] != .data[["SNP_B"]]) %>%
		dplyr::mutate(PHASE=gsub("/", "", .data[["PHASE"]])) %>%
		dplyr::filter(nchar(.data[["PHASE"]]) == 4) %>%
		dplyr::mutate(SNP_A = paste0("rs", .data[["SNP_A"]]), SNP_B = paste0("rs", .data[["SNP_B"]]))

	temp <- do.call(rbind, strsplit(ld[["PHASE"]], "")) %>% dplyr::as_tibble(.name_repair="minimal")
	names(temp) <- c("A1", "B1", "A2", "B2")
	ld <- cbind(ld, temp) %>% dplyr::as_tibble(.name_repair="minimal")
	ld <- dplyr::arrange(ld, dplyr::desc(abs(.data[["R"]])))
	message("Found ", nrow(ld), " proxies")
	RSQLite::dbDisconnect(conn)
	return(ld)
}


#' Extract SNPs from vcf file
#'
#' Finds proxies if necessary
#'
#' @param vcf vcf file name
#' @param rsid list of rs IDs
#' @param bfile ld reference panel (plink)
#' @param proxies ="yes" If SNPs are absent then look for proxies (yes) or not (no). Can also mask all target SNPs and only return proxies (only), for testing purposes
#' @param tag_kb =5000 Proxy parameter
#' @param tag_nsnp =5000 Proxy parameter
#' @param tag_r2 =0.6 Proxy parameter
#' @param threads Number of threads to use (=1)
#' @param rsidx Path to rsidx index 
#' @param dbfile ld tag database (sqlite)
#'
#' @export
#' @return data frame
proxy_match <- function(vcf, rsid, bfile=NULL, proxies="yes", tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1, rsidx=NULL, dbfile=NULL)
{
	if(is.null(dbfile) & is.null(bfile))
	{
		stop('please provide either bfile or dbfile')
	}
	if(!is.null(dbfile) & !is.null(bfile))
	{
		warning("bfile and dbfile both provided. Using dbfile.")
	}
	os <- Sys.info()[['sysname']]
	if(proxies=="yes")
	{
		message("Initial search...")
		o <- query_gwas(vcf, rsid=rsid, rsidx=rsidx)
		missing <- rsid[!rsid %in% names(o)]
		if(length(missing) != 0)
		{
			message("Extracted ", length(rsid) - length(missing), " out of ", length(rsid), " rsids")
			message("Searching for proxies for ", length(missing), " rsids")
			searchspacename <- tempfile()
			if(is.character(vcf))
			{
				if(check_bcftools() & is.null(dbfile))
				{
					message("Determining searchspace...")
					cmd <- paste0(options()[["tools_bcftools"]], " query -f'%ID\n' ", vcf, " > ", searchspacename)
					system(cmd)
					searchspace <- scan(searchspacename, what="character", quiet=TRUE)
				} else {
					searchspace <- NULL
				}				
			} else {
				message("Determining searchspace...")
				searchspace <- names(SummarizedExperiment::rowRanges(vcf))
			}
			message("Proxy lookup...")
			if(is.null(dbfile))
			{
				ld <- get_ld_proxies(missing, bfile, searchspace=searchspace, tag_kb=tag_kb, tag_nsnp=tag_nsnp, tag_r2=tag_r2, threads=threads)
			} else {
				ld <- sqlite_ld_proxies(rsids=missing, dbfile=dbfile, tag_r2=tag_r2)
			}
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
			if(check_bcftools() & is.null(dbfile))
			{
				message("Determining searchspace...")
				cmd <- paste0(options()[["tools_bcftools"]], " query -f'%ID\n' ", vcf, " > ", searchspacename)
				system(cmd)
				searchspace <- scan(searchspacename, what="character")
			} else {
				searchspace <- NULL
			}				
		} else {
			message("Determining searchspace...")
			searchspace <- names(SummarizedExperiment::rowRanges(vcf))
		}
		message("Proxy lookup...")
		if(is.null(dbfile))
		{
			ld <- get_ld_proxies(rsid, bfile, searchspace=searchspace, tag_kb=tag_kb, tag_nsnp=tag_nsnp, tag_r2=tag_r2, threads=threads)			
		} else {
			ld <- sqlite_ld_proxies(rsids=rsid, dbfile=dbfile, tag_r2=tag_r2)
		}
		if(nrow(ld) == 0)
		{
			return(VCF())
		}
	} else if(proxies == "no") {
		o <- query_gwas(vcf, rsid=rsid, rsidx=rsidx)
		return(o)
	} else {
		stop('proxies must be "yes", "no" or "only"')
	}
	if(!is.null(searchspace))
	{
		ld <- ld %>% dplyr::filter(!duplicated(.data[["SNP_A"]]))
	}
	message("Extrating proxies...")
	e <- query_gwas(vcf, rsid=ld[["SNP_B"]], rsidx=rsidx)

	if(is.null(searchspace))
	{
		ld <- subset(ld, ld[["SNP_B"]] %in% names(e)) %>%
			dplyr::filter(!duplicated(.data[["SNP_A"]]))		
	}
	e <- e[names(e) %in% ld[["SNP_B"]], ]
	message("Identified proxies for ", nrow(e), " of ", length(missing), " rsids")
	message("Aligning...")
	index <- match(names(e), ld[["SNP_B"]])
	ld <- ld[index,]
	if(nrow(ld) == 0)
	{
		return(o)
	}
	stopifnot(all(ld[["SNP_B"]] == names(e)))
	sign_index <- GenomicRanges::mcols(SummarizedExperiment::rowRanges(e))[,"REF"] == ld[["B1"]]
	gr <- GenomicRanges::GRanges(ld[["CHR_A"]], IRanges::IRanges(start=ld[["BP_A"]], end=ld[["BP_A"]], names=ld[["SNP_A"]]))
	fixeddat <- S4Vectors::DataFrame(
		REF=Biostrings::DNAStringSet(ld[["A1"]]), 
		ALT=Biostrings::DNAStringSetList(as.list(ld[["A2"]])), 
		QUAL=as.numeric(NA), 
		FILTER="PASS"
	)
	prox <- VariantAnnotation::VCF(
		rowRanges = gr,
		colData = SummarizedExperiment::colData(e),
		fixed = fixeddat,
		info = VariantAnnotation::info(e),
		exptData = list(
			header = VariantAnnotation::header(e)
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

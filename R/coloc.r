#' Perform coloc
#'
#' Perform coloc between two vcfs given a particular locus
#'
#' @param vcf1 VCF object
#' @param vcf2 VCF object
#' @param chrompos chrompos object (see input to parse_chrompos)
#' @param p1 prior prob of assoc in vcf1
#' @param p2 prior prob of assoc in vcf2
#' @param p12 prior prob of joint assoc
#'
#' @export
#' @return Result from coloc
perform_coloc <- function(vcf1, vcf2, p1=1e-4, p2=1e-4, p12=1e-5)
{
	## TODO: binary or quantitative traits
	## TODO: multiallelic variants
	stopifnot(length(vcf1) == length(vcf2))
	tab1 <- vcf1 %>% vcf_to_granges %>% dplyr::as_tibble()
	tab2 <- vcf2 %>% vcf_to_granges %>% dplyr::as_tibble()
	index <- as.character(tab1$REF) == as.character(tab2$REF) &
			as.character(tab1$ALT) == as.character(tab2$ALT) &
			as.character(tab1$seqnames) == as.character(tab2$seqnames) &
			tab1$start == tab2$start

	tab1 <- tab1[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, type = "quant", snp = names(vcf2)[index])}
	tab2 <- tab2[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, type = "quant", snp = names(vcf2)[index])}

	return(coloc::coloc.abf(tab1, tab2, p1=p1, p2=p2, p12=p12))
}



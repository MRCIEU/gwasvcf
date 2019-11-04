#' Perform coloc
#'
#' Perform coloc between two vcfs given a particular locus
#'
#' @param vcf1 VCF object
#' @param vcf2 VCF object
#' @param chrompos chrompos object (see input to parse_chrompos)
#'
#' @export
#' @return Result from coloc
perform_coloc <- function(vcf1, vcf2, p1=1e-4, p2=1e-4, p12=1e-5)
{
	tab1 <- vcf1 %>% vcf_to_granges %>% dplyr::as_tibble()
	tab2 <- vcf2 %>% vcf_to_granges %>% dplyr::as_tibble()
	index <- tab1$REF == tab2$REF &
			tab1$ALT == tab2$ALT &
			tab1$seqnames == tab2$seqnames &
			tab1$start == tab2$start

	tab1 <- tab1[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, type = "quant", snp = names(tab2)[index])}
	tab2 <- tab2[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, type = "quant", snp = names(tab2)[index])}

	return(coloc::coloc.abf(tab1, tab2, p1=p1, p2=p2, p12=p12))
}



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
perform_coloc <- function(vcf1, vcf2, chrompos)
{
	tab1 <- vcf1 %>% vcf_to_granges
	tab2 <- vcf2 %>% vcf_to_granges
	index <- IRanges::mcols(tab1)$REF == IRanges::mcols(tab2)$REF & IRanges::mcols(tab1)$ALT == IRanges::mcols(tab2)$ALT

	tab1 <- tab1[index,] %$% list(pvalues = 10^-LP, N = SS, MAF = AF, beta = ES, varbeta = SE^2, type = "quant", snp = names(tab2)[index])
	tab2 <- tab2[index,] %$% list(pvalues = 10^-LP, N = SS, MAF = AF, beta = ES, varbeta = SE^2, type = "quant", snp = names(tab2)[index])

	return(coloc::coloc.abf(tab1, tab2))
}



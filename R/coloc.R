#' Perform coloc
#'
#' Perform coloc between two vcfs given a particular locus
#'
#' @param vcf1 VCF object
#' @param vcf2 VCF object
#' @param cp =NULL chrompos object (see input to parse_chrompos)
#'
#' @export
#' @return Result from coloc
perform_coloc <- function(vcf1, vcf2, cp=NULL)
{
	if(!is.null(cp))
	{
		vcf1 <- query_gwas(vcf1, cp)
		vcf2 <- query_gwas(vcf2, cp)
	}
	o <- findOverlaps(vcf1, vcf2)
	tab1 <- vcf1[o@from,] %>% vcf_to_granges
	tab2 <- vcf2[o@to,] %>% vcf_to_granges
	index <- mcols(tab1)$REF == mcols(tab2)$REF & mcols(tab1)$ALT == mcols(tab2)$ALT
	tab1 <- tab1[index,] %$% list(pvalues = 10^-LP, N = SS, MAF = AF, beta = ES, varbeta = SE^2, type = "quant", snp = names(tab2)[index])
	tab2 <- tab2[index,] %$% list(pvalues = 10^-LP, N = SS, MAF = AF, beta = ES, varbeta = SE^2, type = "quant", snp = names(tab2)[index])
	return(coloc::coloc.abf(tab1, tab2))
}


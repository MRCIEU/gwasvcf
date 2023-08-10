#  This file contains the functions to create a gwasglue2 SummarySet object



#' Create a SummarySet
#' 
#' Returns a gwasglue2 SummarySet object
#' @param vcf Path or URL to GWAS-VCF file or VCF object e.g. output from VariantAnnotation::readVcf, gwasvcftools::query_vcf or [query_gwas()]
#' @export
gwasvcf_to_summaryset <- function(vcf){
	# get metadata from vcf and create metadata object
	md <- gwasglue2::create_metadata(id = vcf@metadata$header@samples, build = unique(VariantAnnotation::meta(header(vcf))$contig$assembly))
    
	# get summary data and create SummarySet

    s <- vcf %>% 
		vcf_to_tibble() %>% 
		gwasglue2::create_summaryset_from_gwasvcf(metadata = md)

    return(s)
}


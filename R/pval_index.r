#' Create pval index from GWAS-VCF file
#'
#' Create a separate file called <id>.pvali which is used to speed up p-value queries
#'
#' @param vcffile VCF filename
#' @param maximum_pval Maximum p-value to include. Default = 0.05
#' @param indexname index file name to create. Deletes existing file if exists.
#'
#' @export
#' @return NULL
create_pval_index_from_vcf <- function(vcffile, maximum_pval, indexname)
{
    stopifnot(!is.null(options()$tools_bcftools))
    checksqlite3 <- system("which sqlite3")
    if(checksqlite3 != 0) stop("sqlite3 not installed")
	fn <- tempfile()
    id <- VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(vcffile))
    if(length(id) != 1)
    {
        stop("Not implemented for vcf files that don't have a single study")
    }
    cmd <- paste0(options()$tools_bcftools, " query -s ", id, " -i 'FORMAT/LP > ", -log10(maximum_pval), "' -f '%CHROM,%POS,[%LP]\n' ", vcffile, " | sort -r -k 3 > ", fn)
	message("Extracting pval info")
	system(cmd)
    cmd <- c(
        'CREATE TABLE pval_to_coord (chrom TEXT NOT NULL DEFAULT NULL, coord INTEGER NOT NULL DEFAULT NULL, LP REAL NOT NULL DEFAULT 0);',
        '.separator ,',
        paste0('.import ', fn, ' pval_to_coord'),
        'CREATE INDEX idx_LP ON pval_to_coord (LP)'
    )
    print(cmd)
    utils::write.table(cmd, file=paste0(fn, ".sql"), row=F, col=F, qu=F)
	message("Generating index")
	cmd <- paste0("sqlite3 ", indexname, " < ", fn, ".sql")
	unlink(indexname)
	system(cmd)
}

#' Query pval from file using pvali index
#'
#' See create_pvali_index
#'
#' @param pval pval threshold
#' @param vcffile Path to .vcf.gz GWAS summary data file
#' @param id If multiple GWAS datasets in the vcf file, the name (sample ID) from which to perform the filter
#' @param pvali Path to pval index file
#'
#' @export
#' @return vcf object
query_pval_sqlite3 <- function(pval, vcffile, id=NULL, pvali)
{
	out <- query_pvali(pval, pvali)
	return(
		query_gwas(vcffile, chrompos=data.frame(chrom=out$chrom, start=out$coord, end=out$coord), id=id)
	)
}

#' Query pvali
#'
#' @param pval pval threshold
#' @param pvali Path to pval index file
#'
#' @export
#' @return data frame
query_pvali <- function(pval, pvali)
{
	conn <- RSQLite::dbConnect(RSQLite::SQLite(), pvali)
	query <- paste0("SELECT DISTINCT * FROM pval_to_coord WHERE lp >= ", -log10(pval))
	out <- RSQLite::dbGetQuery(conn, query)
	RSQLite::dbDisconnect(conn)
	return(out)
}

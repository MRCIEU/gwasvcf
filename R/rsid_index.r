#' <brief desc>
#'
#' <full description>
#'
#' @param vcf <what param does>
#' @param indexname <what param does>
#'
#' @export
#' @return
create_rsidx_index <- function(vcf, indexname)
{
	fn <- tempfile()
	cmd <- paste0("zcat ", vcf, " | grep -v '#' | awk '{ print substr($3, 3), $1, $2 }' > ", fn, ".txt")
	message("Extracting position info")
	system(cmd)

	cmd <- c(
		'CREATE TABLE rsid_to_coord (rsid INTEGER PRIMARY KEY, chrom TEXT NULL DEFAULT NULL, coord INTEGER NOT NULL DEFAULT 0);',
		'.separator " "',
		paste0('.import ', fn, '.txt rsid_to_coord')
	)
	write.table(cmd, file=paste0(fn, ".sql"))
	cmd <- paste0("sqlite3 ", indexname, " < ", fn, ".sql")
	message("Generating index")
	system(cmd)
}

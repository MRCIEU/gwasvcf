#' <brief desc>
#'
#' <full description>
#'
#' @param vcf VCF filename
#' @param indexname index file name to create. Deletes existing file if exists.
#'
#' @export
#' @return NULL
create_rsidx_index_from_vcf <- function(vcf, indexname)
{
	fn <- tempfile()
	cmd <- paste0("gunzip -c ", vcf, " | grep -v '#' | awk '{ print substr($3, 3), $1, $2 }' > ", fn, ".txt")
	message("Extracting position info")
	system(cmd)

	cmd <- c(
		'CREATE TABLE rsid_to_coord (rsid INTEGER PRIMARY KEY, chrom TEXT NULL DEFAULT NULL, coord INTEGER NOT NULL DEFAULT 0);',
		'.separator " "',
		paste0('.import ', fn, '.txt rsid_to_coord')
	)
	write.table(cmd, file=paste0(fn, ".sql"), row=F, col=F, qu=F)
	message("Generating index")
	cmd <- paste0("sqlite3 ", indexname, " < ", fn, ".sql")
	unlink(indexname)
	system(cmd)
}

#' Create new index from existing index using a subset of rsids
#'
#' @param rsid Vector of rsids
#' @param rsidx Existing index
#' @param newindex New index (Note: will delete existing file if exists)
#'
#' @export
#' @return NULL, creates new index file
create_rsidx_sub_index <- function(rsid, rsidx, newindex)
{
	out <- query_rsidx(rsid, rsidx)
	unlink(newindex)
	conn <- RSQLite::dbConnect(RSQLite::SQLite(), newindex)
	RSQLite::dbWriteTable(conn, "rsid_to_coord", out)
	RSQLite::dbExecute(conn, "CREATE INDEX rsid on rsid_to_coord(rsid);")
	RSQLite::dbDisconnect(conn)
}

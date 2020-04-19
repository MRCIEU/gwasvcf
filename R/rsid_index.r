#' <brief desc>
#'
#' <full description>
#'
#' @param vcf <what param does>
#' @param indexname <what param does>
#'
#' @export
#' @return
create_rsidx_index_from_vcf <- function(vcf, indexname)
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
	write.table(cmd, file=paste0(fn, ".sql"), row=F, col=F, qu=F)
	cmd <- paste0("sqlite3 ", indexname, " < ", fn, ".sql")
	message("Generating index")
	system(cmd)
}

#' @export
create_rsidx_sub_index <- function(rsid, rsidx, newindex)
{
	conn <- RSQLite::dbConnect(RSQLite::SQLite(), rsidx)
	numid <- gsub("rs", "", rsid) %>% paste(., collapse=",")
	query <- paste0("SELECT DISTINCT * FROM rsid_to_coord WHERE rsid IN (", numid, ")")
	out <- RSQLite::dbGetQuery(conn, query)
	dbDisconnect(conn)
	unlink(newindex)
	conn <- dbConnect(RSQLite::SQLite(), newindex)
	dbWriteTable(conn, "rsid_to_coord", out)
	dbExecute(conn, "CREATE INDEX rsid on rsid_to_coord(rsid);")
	dbDisconnect(conn)
}
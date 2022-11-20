#' Create RSID index from VCF
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
	utils::write.table(cmd, file=paste0(fn, ".sql"), row=F, col=F, qu=F)
	message("Generating index")
	cmd <- paste0("sqlite3 ", indexname, " < ", fn, ".sql")
	unlink(indexname)
	system(cmd)
}

#' Create new index from existing index using a subset of rsids
#'
#' Note this requires a modified version of plink that allows ld-window-r2 flag for --r option.
#' Available here: https://github.com/explodecomputer/plink-ng
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



#' Create LD reference sqlite database for tags
#'
#' This is used for looking up proxies
#'
#' @param bfile path to plink file
#' @param dbname dbname to produce (overwrites existing if exists)
#' @param tag_r2 minimum tag r2
#'
#' @export
#' @return NULL
create_ldref_sqlite <- function(bfile, dbname, tag_r2=0.6)
{
	stopifnot(check_plink())
	message("identifying indels to remove")
	cmd <- paste0("awk '{ if (length($5) != 1 || length($6) != 1) { print $2 }}' ", bfile, ".bim > ", bfile, ".indels")
	system(cmd)

	message("calculating ld tags")
	cmd <- paste0(options()[["tools_plink"]], " --bfile ", bfile, " --keep-allele-order --exclude ", bfile, ".indels --r in-phase with-freqs gz --out ", bfile, " --ld-window-kb 250 --ld-window 1000 --ld-window-r2 ",  tag_r2)
	system(cmd)

	message("formatting")
	cmd <- paste0("gunzip -c ", bfile, ".ld.gz | awk 'BEGIN {OFS=\",\"}  { if(NR != 1) { print substr($3, 3), $1, $2, $4, substr($7, 3), $5, $6, $9, $8, $10 }}' > ", bfile, ".ld.tab")
	system(cmd)

	message("creating sqlite db")
	cmd <- c(
		'CREATE TABLE tags (',
		'	SNP_A INTEGER NOT NULL, ',
		'	CHR_A TEXT NULL DEFAULT NULL, ',
		'	BP_A INTEGER NOT NULL,',
		'	MAF_A REAL NOT NULL,',
		'	SNP_B INTEGER NOT NULL, ',
		'	CHR_B TEXT NULL DEFAULT NULL, ',
		'	BP_B INTEGER NOT NULL,',
		'	MAF_B REAL NOT NULL,',
		'	PHASE TEXT NOT NULL,',
		'	R REAL NOT NULL',
		');',
		'CREATE INDEX SNP_A_INDEX ON tags(SNP_A);',
		'.separator ","',
		paste0(".import ", bfile, ".ld.tab tags")
	)
	unlink(paste0(bfile, ".ld.sqlite"))
	utils::write.table(cmd, file=paste0(bfile, ".ld.sqlite"), row=F, col=F, qu=F)
	unlink(dbname)
	cmd <- paste0("sqlite3 ", dbname, " < ", bfile, ".ld.sqlite")
	system(cmd)
	unlink(paste0(bfile, ".ld.tab"))
	unlink(paste0(bfile, ".ld.gz"))
	unlink(paste0(bfile, ".ld.sqlite"))
	# unlink(paste0(bfile, ".indels"))
}


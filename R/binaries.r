 #' Check if the tools_bcftools option is set
#'
#' See set_bcftools() for more information
#'
#'
#' @export
#' @return TRUE or FALSE
check_bcftools <- function()
{
	if(is.null(options()$tools_bcftools))
	{
		message("'tools_bcftools' option is not set, using native read which may be substantially slower. See 'set_bcftools' for information.")
		return(FALSE)
	} else if(!file.exists(options()$tools_bcftools)) {
		message("'tools_bcftools' option does not point to an existing file, using native read which may be substantially slower. See 'set_bcftools' for information.")
		return(FALSE)
	} else {
		return(TRUE)
	}
}


#' Check if the tools_plink option is set
#'
#' See set_plink() for more information
#'
#'
#' @export
#' @return TRUE or FALSE
check_plink <- function()
{
	if(is.null(options()$tools_plink))
	{
		message("'tools_plink' option is not set. See 'set_plink' for information.")
		return(FALSE)
	} else if(!file.exists(options()$tools_plink)) {
		message("'tools_plink' option does not point to an existing file. See 'set_plink' for information.")
		return(FALSE)
	} else {
		return(TRUE)
	}
}

#' Set bcftools binary location
#'
#'
#' @param path If "" (default), then will use the explodecomputer/genetics.binaRies to get binaries that are appropriate for the detected operating system. Otherwise, provide the path to the bcftools binary. If NULL then will set the option to NULL.
#'
#' @export
#' @return NULL, sets option 'tools_bcftools'
set_bcftools <- function(path="")
{
	if(is.null(path))
	{
		options(tools_bcftools = NULL)
	} else if(path == "")
	{
		a <- require(genetics.binaRies)
		if(a)
		{
			message("Path not provided, using binaries in the explodecomputer/genetics.binaRies package")
			options(tools_bcftools = genetics.binaRies::get_bcftools_binary())
		} else {
			stop("Please provide a path to bcftools binary or run devtools::install_github('explodecomputer/genetics.binaRies')")
		}
	} else {
		options(tools_bcftools = path)
	}
}

#' Set plink binary location
#'
#'
#' @param path If "" (default), then will use the explodecomputer/genetics.binaRies to get binaries that are appropriate for the detected operating system. Otherwise, provide the path to the plink binary. If NULL then will set the option to NULL.
#'
#' @export
#' @return NULL, sets option 'tools_plink'
set_plink <- function(path="")
{
	if(is.null(path))
	{
		options(tools_plink = NULL)
	} else if(path == "")
	{
		a <- require(genetics.binaRies)
		if(a)
		{
			message("Path not provided, using binaries in the explodecomputer/genetics.binaRies package")
			options(tools_plink = genetics.binaRies::get_plink_binary())
		} else {
			stop("Please provide a path to plink binary or run devtools::install_github('explodecomputer/genetics.binaRies')")
		}
	} else {
		options(tools_plink = path)
	}
}

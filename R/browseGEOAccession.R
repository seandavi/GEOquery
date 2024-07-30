#' The URL for a GEO accession
#'
#' Sometimes, you just need the URL for a GEO accession. This function
#' returns the URL for a given GEO accession number that
#' can be used to access the GEO page for that accession.
#'
#' @param geo A GEO accession number
#'
#' @return A character vector with the URL for the GEO accession
#'
#' @seealso \code{\link{browseGEOAccession}}
#'
#' @examples
#'
#' urlForAccession("GSE262484")
#'
#' @export
urlForAccession <- function(geo) {
  paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", geo)
}

#' Open the GEO page for a given accession
#'
#' Sometimes, you just need to see the GEO website page for a
#' GEO accession. This function opens the GEO page for a given
#' accession number in the default browser.
#'
#' @param geo A GEO accession number
#'
#' @seealso \code{\link{urlForAccession}}
#'
#' @examples
#' \dontrun{
#' browseGEOAccession("GSE262484")
#' }
#'
#' @export
browseGEOAccession <- function(geo) {
  browseURL(urlForAccession(geo))
}

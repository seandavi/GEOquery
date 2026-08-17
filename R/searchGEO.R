#' Search GEO database
#'
#' This function searches the [GDS](https://www.ncbi.nlm.nih.gov/gds)
#' database, and return a data.frame for all the search results.
#'
#' The NCBI allows users to access more records (10 per second) if they register
#' for and use an API key. [set_entrez_key][rentrez::set_entrez_key] function
#' allows users to set this key for all calls to rentrez functions during a
#' particular R session. You can also set an environment variable `ENTREZ_KEY`
#' by [Sys.setenv][base::Sys.setenv].  Once this value is set to your key
#' rentrez will use it for all requests to the NCBI. Details see
#' <https://docs.ropensci.org/rentrez/articles/rentrez_tutorial.html#rate-limiting-and-api-keys>
#'
#' @param query character, the search term. The NCBI uses a search term syntax
#' which can be associated with a specific search field with square brackets.
#' So, for instance "Homo sapiens\[ORGN\]" denotes a search for `Homo sapiens`
#' in the “Organism” field. Details see
#' <https://www.ncbi.nlm.nih.gov/geo/info/qqtutorial.html>. The names and
#' definitions of these fields can be identified using
#' [searchFieldsGEO].
#'
#' @seealso [searchFieldsGEO]
#'
#' @param step the number of records to fetch from the database each time. You
#' may choose a smaller value if failed.
#'
#' @return a data.frame with one row per matching GEO record and columns
#' `Accession`, `Title`, `Summary`, `Organism`, `Type`, `GPL`, `n_samples`,
#' `PDAT`, `suppFile`, `FTPLink`, `SeriesTitle`, `entryType`, and `ID` (the
#' Entrez UID). An empty query returns a zero-row data.frame with those columns.
#'
#' @examples
#' \dontrun{
#' searchGEO("diabetes[ALL] AND Homo sapiens[ORGN] AND GSE[ETYP]")
#' }
#'
#' @importFrom rentrez entrez_search
#' @export
searchGEO <- function(query, step = 500L) {
  search_res <- rentrez::entrez_search(
    "gds", query,
    use_history = TRUE, retmax = 0L
  )
  count <- search_res$count
  if (!count) {
    return(.empty_gds_summary())
  }
  # esummary retstart is 0-based.
  seq_starts <- seq(0L, count - 1L, by = step)
  records <- vector("list", length(seq_starts))
  for (i in seq_along(seq_starts)) {
    json <- .fetch_gds_esummary_json(
      web_history = search_res$web_history,
      retstart = seq_starts[[i]], retmax = step
    )
    records[[i]] <- .parse_gds_esummary_json(json)
    Sys.sleep(1L)
  }
  as.data.frame(
    data.table::rbindlist(records, use.names = TRUE, fill = TRUE)
  )
}


#' Provide a list of possible search fields for GEO search
#'
#' @returns a data.frame with names of possible search fields for GEO search
#' as well as descriptions, data types, etc. for each field. Fields are
#' in rows and their properties are in columns.
#'
#' @seealso \code{\link{searchGEO}}
#'
#' @examples
#' searchFieldsGEO()
#'
#' @importFrom rentrez entrez_db_searchable
#' @export
searchFieldsGEO <- function() {
  res <- do.call(
    rbind,
    rentrez::entrez_db_searchable("gds")
  ) |> data.frame()
  rownames(res) <- NULL
  res
}


# Fetch one batch of GDS esummary records as JSON (esummary version 2.0). Uses
# the Entrez history from a prior entrez_search(use_history=TRUE). Kept separate
# from the parser so the parser is testable offline; this is the only piece that
# touches the network. An ENTREZ_KEY, if set, raises the NCBI rate limit.
.fetch_gds_esummary_json <- function(web_history, retstart, retmax,
    timeout = getOption("GEOquery.download.timeout", 300)) {
  req <- .geo_request(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", timeout
  ) |>
    httr2::req_url_query(
      db = "gds", version = "2.0", retmode = "json",
      WebEnv = web_history$WebEnv, query_key = web_history$QueryKey,
      retstart = retstart, retmax = retmax
    )
  key <- Sys.getenv("ENTREZ_KEY")
  if (nzchar(key)) {
    req <- httr2::req_url_query(req, api_key = key)
  }
  httr2::resp_body_string(httr2::req_perform(req))
}

# Parse a GDS esummary (version 2.0) JSON response into a tidy data.frame, one
# row per record. Offline-testable core of searchGEO().
#' @importFrom jsonlite fromJSON
.parse_gds_esummary_json <- function(json) {
  result <- jsonlite::fromJSON(json, simplifyVector = FALSE)$result
  uids <- unlist(result$uids)
  if (is.null(result) || !length(uids)) {
    return(.empty_gds_summary())
  }
  rows <- lapply(uids, function(u) .gds_summary_row(result[[u]]))
  as.data.frame(
    data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  )
}

# Map one esummary record (a named list) to a one-row data.frame with stable,
# human-meaningful columns. Missing fields become NA.
.gds_summary_row <- function(rec) {
  g <- function(field) {
    v <- rec[[field]]
    if (is.null(v) || !length(v)) NA_character_ else as.character(v)[[1L]]
  }
  data.frame(
    Accession = g("accession"),
    Title = g("title"),
    Summary = g("summary"),
    Organism = g("taxon"),
    Type = g("gdstype"),
    GPL = g("gpl"),
    n_samples = suppressWarnings(as.integer(g("n_samples"))),
    PDAT = g("pdat"),
    suppFile = g("suppfile"),
    FTPLink = g("ftplink"),
    SeriesTitle = g("seriestitle"),
    entryType = g("entrytype"),
    ID = g("uid"),
    stringsAsFactors = FALSE
  )
}

.empty_gds_summary <- function() {
  .gds_summary_row(list())[0L, , drop = FALSE]
}

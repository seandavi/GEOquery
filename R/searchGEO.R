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
#' @return a data.frame contains the search results
#'
#' @examples
#' \dontrun{
#' searchGEO("diabetes[ALL] AND Homo sapiens[ORGN] AND GSE[ETYP]")
#' }
#'
#' @export
searchGEO <- function(query, step = 500L) {
  records_num <- rentrez::entrez_search(
    "gds", query,
    retmax = 0L
  )$count
  seq_starts <- seq(1L, records_num, step)
  records <- character(length(seq_starts))
  search_res <- rentrez::entrez_search(
    "gds", query,
    use_history = TRUE, retmax = 0L
  )
  for (i in seq_along(seq_starts)) {
    records[[i]] <- rentrez::entrez_fetch(
      db = "gds", web_history = search_res$web_history,
      rettype = "summary", retmode = "text",
      retmax = step, retstart = seq_starts[[i]]
    )
    Sys.sleep(1L)
  }
  records <- strsplit(
    gsub("^\\n|\\n$", "", paste0(records, collapse = "")),
    "\\n\\n"
  )[[1L]]
  name_value_pairs <- parse_name_value_pairs(preprocess_records(records))
  data.table::setDF(name_value_pairs)
  name_value_pairs
}


#' Provide a list of possible search fields for GEO search
#'
#' @import rentrez
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
#' @export
searchFieldsGEO <- function() {
  res <- do.call(
    rbind,
    rentrez::entrez_db_searchable("gds")
  ) |> data.frame()
  rownames(res) <- NULL
  res
}



# this function just processed GEO searched results returned by `entrez_fetch`
# into key-values paris
preprocess_records <- function(x) {
  x <- sub("^\\d+\\.", "Title:", x, perl = TRUE)
  x <- sub(
    "\\n\\(Submitter supplied\\)\\s*",
    "\nSummary: ", x,
    perl = TRUE
  )
  x <- gsub(
    "\\s*\\n?(Platform|Dataset)s?\\s*:\\s*",
    "\n\\1s: ", x,
    perl = TRUE
  )
  x <- sub("\\tID:\\s*", "\nID: ", x, perl = TRUE)
  x <- sub(
    "\\n?\\s*((?:\\s*\\d+(?:\\s*related)?\\s*(?:DataSet|Platform|Sample|Serie)s?)+)([^:])",
    "\nContains: \\1\\2", x,
    perl = TRUE
  )
  x <- gsub(":\\t+", ": ", x, perl = TRUE)
  x <- gsub("\\t\\t+", " ", x, perl = TRUE)
  strsplit(x, "\\n", perl = TRUE)
}

# parse key-value pairs separeted by ":". For a list of key-value pairs
# characters (like: `list(c("a:1", "b:2"), c("a:3", "b:4"))`), this function
# simply cleans those up and transforms the list into a list object, the names
# of returned value is the unique keys in the pairs, the element of the returned
# list is the values in the paris.
# See parse_name_value_pairs(list(c("a:1", "b:2"), c("a:3", "b:4")))
#' @return a list, every element of which corresponds to each key-value pairs
#' group by key in the paris.
#' @noRd
parse_name_value_pairs <- function(chr_list, sep = ":") {
  .characteristic_list <- lapply(chr_list, function(x) {
    if (!length(x)) {
      return(data.table::data.table())
    }
    # Don't use `data.table::tstrsplit`, as it will split string into three
    # or more elements.
    name_value_pairs <- data.table::transpose(
      str_split(x, paste0("(\\s*+)", sep, "(\\s*+)"))
    )
    res <- as.list(name_value_pairs[[2L]])
    names(res) <- name_value_pairs[[1L]]
    data.table::setDT(res)
    res
  })
  characteristic_dt <- data.table::rbindlist(
    .characteristic_list,
    use.names = TRUE, fill = TRUE
  )
  data.table::setnames(characteristic_dt, make.unique)

  # parse text into corresponding atomic vector mode
  lapply(characteristic_dt, function(x) {
    data.table::fread(
      text = x, sep = "", header = FALSE,
      strip.white = TRUE, blank.lines.skip = FALSE, fill = TRUE
    )[[1L]]
  })
}

# split string based on pattern, Only split once, Return a list of character,
# the length of every element is two
str_split <- function(string, pattern, ignore.case = FALSE) {
  regmatches(
    string,
    regexpr(pattern, string,
      perl = TRUE, fixed = FALSE,
      ignore.case = ignore.case
    ),
    invert = TRUE
  )
}

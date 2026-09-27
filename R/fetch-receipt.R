# Optional retrieval receipts (#221; ADR-0009). The active collector belongs to
# one outer getGEO()/getGEOSuppFiles() call, including its nested GPL downloads.
# It is cleared on exit, including errors. No work is done when disabled.
.geo_receipt_state <- new.env(parent = emptyenv())

.receipt_error <- function(message) {
  .geoquery_abort(message, class = "geoquery_receipt_error")
}

.receipt_path_exists <- function(path) {
  if (file.exists(path)) {
    return(TRUE)
  }
  # file.exists() misses dangling symlinks; preserve the link itself.
  link <- Sys.readlink(path)
  !is.na(link) && nzchar(link)
}

.receipt_start <- function(operation) {
  if (!is.null(.geo_receipt_state$active)) {
    return(NULL)
  }
  path <- getOption("GEOquery.fetch.receipt")
  if (is.null(path) || identical(path, FALSE)) {
    return(NULL)
  }
  if (!is.character(path) || length(path) != 1L ||
        is.na(path) || !nzchar(path)) {
    .receipt_error(
      "'GEOquery.fetch.receipt' must be NULL, FALSE or a single output path."
    )
  }
  path <- path.expand(path)
  parent <- dirname(path)
  if (!dir.exists(parent) || file.access(parent, 2L) != 0L) {
    .receipt_error("The receipt's parent directory must exist and be writable.")
  }
  path <- file.path(
    normalizePath(parent, winslash = "/", mustWork = TRUE), basename(path)
  )
  if (.receipt_path_exists(path)) {
    .receipt_error("The receipt path already exists; choose a new path.")
  }
  receipt <- new.env(parent = emptyenv())
  receipt$path <- path
  receipt$operation <- operation
  receipt$resources <- list()
  .geo_receipt_state$active <- receipt
  receipt
}

.receipt_stop <- function(receipt) {
  if (!is.null(receipt) && identical(.geo_receipt_state$active, receipt)) {
    .geo_receipt_state$active <- NULL
  }
  invisible(NULL)
}

# Only GEO routing parameters are retained. In particular, reviewer tokens,
# arbitrary signed-query credentials, userinfo and fragments never reach disk.
.receipt_url <- function(url) {
  if (is.null(url)) {
    return(NULL)
  }
  parsed <- httr2::url_parse(url)
  parsed$username <- NULL
  parsed$password <- NULL
  parsed$fragment <- NULL
  parsed$query <- parsed$query[intersect(
    names(parsed$query), c("acc", "targ", "form", "view", "tool")
  )]
  httr2::url_build(parsed)
}

.receipt_accession <- function(url) {
  if (is.null(url)) {
    return(NULL)
  }
  # FTP directory stubs such as GSE123nnn are not accessions.
  matches <- regmatches(url, regexpr(
    "(?<![[:alnum:]])G(SE|SM|PL|DS)[0-9]+(?![[:alnum:]])", url, perl = TRUE
  ))
  if (length(matches)) matches else NULL
}

.receipt_record <- function(path, url = NULL, response = NULL,
                            source = "local") {
  receipt <- .geo_receipt_state$active
  if (is.null(receipt)) {
    return(invisible(NULL))
  }
  key <- normalizePath(path, winslash = "/", mustWork = TRUE)
  # The same GPL can be reused by several matrices. Preserve the original
  # download observation rather than replace it with a later local reuse.
  if (!is.null(receipt$resources[[key]])) {
    return(invisible(NULL))
  }
  url <- .receipt_url(url)
  observed_url <- if (source == "local") NULL else url
  endpoint <- observed_url
  if (!is.null(endpoint)) {
    parsed <- httr2::url_parse(endpoint)
    parsed$query <- NULL
    endpoint <- httr2::url_build(parsed)
  }
  receipt$resources[[key]] <- list(
    accession = .receipt_accession(url),
    url = observed_url,
    endpoint = endpoint,
    expected_url = if (source == "local") url else NULL,
    file = basename(path),
    source = source,
    bytes = unname(file.info(path)$size),
    md5 = .file_md5(path),
    http_last_modified = if (is.null(response)) NULL else
      httr2::resp_header(response, "last-modified"),
    geo_last_update = NULL,
    geoquery_version = as.character(utils::packageVersion("GEOquery")),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
  invisible(NULL)
}

# Read update dates from the object already parsed by GEOquery. Supplementary
# formats and old HTTP/cache metadata remain unknown; no extra requests.
.receipt_metadata <- function(path, object) {
  receipt <- .geo_receipt_state$active
  if (is.null(receipt)) {
    return(invisible(NULL))
  }
  key <- normalizePath(path, winslash = "/", mustWork = TRUE)
  record <- receipt$resources[[key]]
  if (is.null(record)) {
    return(invisible(NULL))
  }
  metadata <- if (methods::is(object, "GEOData")) {
    Meta(object)
  } else if (methods::is(object, "ExpressionSet")) {
    Biobase::experimentData(object)@other
  } else {
    return(invisible(NULL))
  }
  scalar <- function(x) {
    if (is.null(x) || !length(x) || is.na(x[[1L]])) {
      NULL
    } else {
      as.character(x[[1L]])
    }
  }
  record["geo_last_update"] <- list(scalar(metadata$last_update_date))
  if (is.null(record$accession)) {
    record["accession"] <- list(scalar(metadata$geo_accession))
  }
  receipt$resources[[key]] <- record
  invisible(NULL)
}

.receipt_write <- function(receipt) {
  if (is.null(receipt)) {
    return(invisible(NULL))
  }
  # Refuse overwrites, including collisions with a file just downloaded by
  # this call. A failed call never creates a completed receipt.
  if (.receipt_path_exists(receipt$path)) {
    .receipt_error("The receipt path already exists; choose a new path.")
  }
  temporary <- tempfile(".geoquery-receipt-", tmpdir = dirname(receipt$path))
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(
    list(
      schema_version = 1L,
      operation = receipt$operation,
      resources = unname(receipt$resources)
    ),
    temporary, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"
  )
  if (!file.rename(temporary, receipt$path)) {
    .receipt_error("Could not move the completed receipt to its output path.")
  }
  invisible(NULL)
}

# Structured error conditions for GEOquery.
#
# Every GEOquery error inherits from `geoquery_error` (itself an rlang error),
# so downstream code can `tryCatch()` on the parent class or on a specific
# subclass. Each constructor may carry structured data (accession, url, status,
# fname) and an optional `parent` condition to preserve the original cause.
#
# These are internal helpers; the documented contract is the set of condition
# classes:
#   geoquery_error              -- parent of all GEOquery conditions
#     geoquery_bad_accession    -- malformed/unrecognized accession (carries $accession)
#     geoquery_private_accession-- non-public/HTML response (carries $accession)
#     geoquery_download_error   -- a download failed (carries $url, $status)
#     geoquery_parse_error      -- parsing a GEO file failed (carries $fname)

#' @importFrom rlang abort
.geoquery_abort <- function(message, class, ..., parent = NULL) {
    rlang::abort(message, class = c(class, "geoquery_error"), ..., parent = parent)
}

.abort_bad_accession <- function(accession) {
    .geoquery_abort(
        sprintf(
            "'%s' is not a recognized GEO accession (expected GSE, GSM, GPL, or GDS followed by digits).",
            accession
        ),
        class = "geoquery_bad_accession", accession = accession
    )
}

.abort_private_accession <- function(accession = NULL) {
    .geoquery_abort(
        paste0(
            "The downloaded content looks like an HTML page, not GEO data. ",
            "The accession may be private, embargoed, not yet public, or may not exist."
        ),
        class = "geoquery_private_accession", accession = accession
    )
}

.abort_download <- function(message, url = NULL, status = NULL, parent = NULL) {
    .geoquery_abort(
        message,
        class = "geoquery_download_error", url = url, status = status, parent = parent
    )
}

.abort_parse <- function(message, fname = NULL, parent = NULL) {
    .geoquery_abort(
        message,
        class = "geoquery_parse_error", fname = fname, parent = parent
    )
}

# Reachability helper for tests that hit NCBI GEO specifically. General
# internet access does not guarantee the GEO hosts are up, so probe them
# directly. Mirrors the BiocPkgTools `skip_if_bioc_offline()` pattern.
#
# This complements `skip_if_no_integration()` (helper-integration.R), which
# gates the live-network suite on GEOQUERY_INTEGRATION for CI determinism.
# Use `skip_if_geo_offline()` for tests that should run whenever GEO is
# reachable, independent of the CI env flag.
.geo_reachable <- function(timeout_sec = 5) {
    if (!curl::has_internet()) {
        return(FALSE)
    }
    probe <- function(url) {
        tryCatch(
            {
                req <- httr2::request(url)
                req <- httr2::req_options(req, timeout_ms = timeout_sec * 1000)
                req <- httr2::req_method(req, "HEAD")
                req <- httr2::req_error(req, is_error = function(resp) FALSE)
                resp <- httr2::req_perform(req)
                httr2::resp_status(resp) < 400
            },
            error = function(e) FALSE
        )
    }
    # Either GEO host being up is enough to exercise the download paths.
    isTRUE(probe("https://www.ncbi.nlm.nih.gov/geo/")) ||
        isTRUE(probe("https://ftp.ncbi.nlm.nih.gov/geo/"))
}

#' Skip a test unless NCBI GEO is reachable
#'
#' Use for tests that hit GEO-hosted resources when you want them to run
#' whenever the network permits (e.g. local runs), rather than only under
#' the CI integration flag.
skip_if_geo_offline <- function() {
    testthat::skip_if_not(
        .geo_reachable(), "NCBI GEO is not reachable"
    )
}

# Live-network ("integration") tests hit NCBI GEO and are non-deterministic in
# CI. They run only when GEOQUERY_INTEGRATION=true (set by the scheduled
# integration workflow), and never on CRAN. Offline/synthetic tests are
# unaffected. See issue #169.
skip_if_no_integration <- function() {
    testthat::skip_on_cran()
    if (!identical(Sys.getenv("GEOQUERY_INTEGRATION"), "true")) {
        testthat::skip("live-network integration test; set GEOQUERY_INTEGRATION=true to run")
    }
}

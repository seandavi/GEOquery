# Offline tests for private-record reviewer-token support (#154, ADR-0007).
# Private GEO records can only be reached via acc.cgi with a `token=` query
# parameter, so these assert URL construction (no network) by mocking the
# download layer to capture the URL that would be fetched.

test_that(".append_token appends a token only when supplied (#154)", {
    expect_equal(GEOquery:::.append_token("https://x/acc.cgi?acc=GSE1"), "https://x/acc.cgi?acc=GSE1")
    expect_null(attr(GEOquery:::.append_token("https://x/acc.cgi?acc=GSE1", NULL), "token"))
    # URL already has a query string -> use '&'
    expect_equal(
        GEOquery:::.append_token("https://x/acc.cgi?acc=GSE1", "abc"),
        "https://x/acc.cgi?acc=GSE1&token=abc"
    )
    # URL without a query string -> use '?'
    expect_equal(
        GEOquery:::.append_token("https://x/plain", "abc"),
        "https://x/plain?token=abc"
    )
})

# Capture the URL getGEOfile() would download by mocking the internal
# downloadFile() to record its first argument and create the destfile.
capture_geofile_url <- function(GEO, ..., token = NULL) {
    seen <- new.env()
    testthat::local_mocked_bindings(
        downloadFile = function(url, destfile, ...) {
            seen$url <- url
            file.create(destfile)
            invisible(TRUE)
        },
        .package = "GEOquery"
    )
    dest <- tempfile()
    dir.create(dest)
    GEOquery:::getGEOfile(GEO, destdir = dest, token = token, ...)
    seen$url
}

test_that("a token routes a GSE through acc.cgi, not the FTP tree (#154)", {
    url <- capture_geofile_url("GSE123456", token = "SECRET")
    expect_match(url, "acc\\.cgi", fixed = FALSE)
    expect_match(url, "token=SECRET", fixed = TRUE)
    expect_match(url, "acc=GSE123456", fixed = TRUE)
    expect_false(grepl("ftp\\.ncbi", url))
})

test_that("without a token a GSE still uses the FTP family SOFT (#154 regression)", {
    url <- capture_geofile_url("GSE123456", token = NULL)
    expect_match(url, "ftp\\.ncbi.*_family\\.soft\\.gz", fixed = FALSE)
    expect_false(grepl("token=", url))
})

test_that("a token is appended for GSM and GPL acc.cgi requests (#154)", {
    gsm <- capture_geofile_url("GSM999", token = "T0K")
    expect_match(gsm, "acc\\.cgi", fixed = FALSE)
    expect_match(gsm, "token=T0K", fixed = TRUE)

    gpl <- capture_geofile_url("GPL570", token = "T0K")
    expect_match(gpl, "acc\\.cgi", fixed = FALSE)
    expect_match(gpl, "token=T0K", fixed = TRUE)
})

test_that("getGEO() validates the token argument (#154)", {
    expect_error(getGEO("GSE1", token = c("a", "b")), "single character string")
    expect_error(getGEO("GSE1", token = 42), "single character string")
})

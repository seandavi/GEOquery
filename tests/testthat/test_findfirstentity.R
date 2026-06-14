# Offline tests for findFirstEntity() hardening (#58). No network: feed crafted
# content through a textConnection.

ff <- function(txt) {
    con <- textConnection(txt)
    on.exit(close(con))
    GEOquery:::findFirstEntity(con)
}

test_that("findFirstEntity errors clearly on an HTML / non-public page (#58)", {
    html <- c(
        "<!DOCTYPE html>", "<html><body>",
        "This record could not be found or is not currently public.",
        "</body></html>"
    )
    expect_error(ff(html), "private|HTML|public|not exist")
})

test_that("findFirstEntity returns 0 for a series-matrix file", {
    sm <- c(
        "!Series_title\tMy series",
        "!Series_geo_accession\tGSE1",
        "!series_matrix_table_begin"
    )
    expect_equal(ff(sm), 0)
})

test_that("findFirstEntity detects a SOFT sample entity", {
    soft <- c("^SAMPLE = GSM12345", "!Sample_title = x")
    res <- ff(soft)
    expect_equal(res[1], "sample")
    expect_equal(res[2], "GSM12345")
})

test_that("findFirstEntity tolerates multiple !Series_title lines (#58)", {
    # Two matching lines in one chunk previously errored on `if (!grepl(...))`
    # because the condition had length > 1 (R >= 4.2).
    expect_equal(ff(c("!Series_title\tA", "!Series_title\tB")), 0)
})

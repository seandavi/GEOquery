# Offline tests for the structured error conditions (#170).

test_that("condition constructors set the class hierarchy and carry data", {
    e <- tryCatch(GEOquery:::.abort_bad_accession("FOO123"), error = function(e) e)
    expect_s3_class(e, "geoquery_bad_accession")
    expect_s3_class(e, "geoquery_error")
    expect_equal(e$accession, "FOO123")

    d <- tryCatch(
        GEOquery:::.abort_download("boom", url = "http://x", status = 404L),
        error = function(e) e
    )
    expect_s3_class(d, "geoquery_download_error")
    expect_s3_class(d, "geoquery_error")
    expect_equal(d$status, 404L)
    expect_equal(d$url, "http://x")

    p <- tryCatch(GEOquery:::.abort_parse("bad", fname = "f.soft"), error = function(e) e)
    expect_s3_class(p, "geoquery_parse_error")
    expect_equal(p$fname, "f.soft")
})

test_that("findFirstEntity raises a typed private-accession condition (#58, #170)", {
    con <- textConnection(c("<!DOCTYPE html>", "This record is not currently public"))
    on.exit(close(con))
    expect_error(
        GEOquery:::findFirstEntity(con),
        class = "geoquery_private_accession"
    )
})

# Offline tests for url_join() (#131): joining URLs must not collapse the
# scheme's '//' (as file.path/fs::path do) nor produce a double slash.

test_that("url_join keeps a single slash and preserves the scheme (#131)", {
    base <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48350/suppl/"
    expect_equal(
        GEOquery:::url_join(base, "GSE48350_RAW.tar"),
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48350/suppl/GSE48350_RAW.tar"
    )
    # scheme '//' intact, no '//' after it
    expect_false(grepl("[^:]//", GEOquery:::url_join(base, "x")))
})

test_that("url_join works whether or not the base has a trailing slash", {
    expect_equal(GEOquery:::url_join("https://a.b/c/", "d"), "https://a.b/c/d")
    expect_equal(GEOquery:::url_join("https://a.b/c", "d"), "https://a.b/c/d")
})

test_that("url_join strips a leading slash on the path and is vectorized", {
    expect_equal(GEOquery:::url_join("https://a.b/c/", "/d"), "https://a.b/c/d")
    expect_equal(
        GEOquery:::url_join("https://a.b/c/", c("d", "e")),
        c("https://a.b/c/d", "https://a.b/c/e")
    )
})

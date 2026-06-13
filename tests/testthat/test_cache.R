# Offline tests for the BiocFileCache-backed download cache (#171). All cache
# operations are local (no network).

test_that("the download cache stores and retrieves by URL (#171)", {
    skip_if_not_installed("BiocFileCache")
    old <- getOption("GEOquery.cache.path")
    on.exit(options(GEOquery.cache.path = old), add = TRUE)
    options(GEOquery.cache.path = tempfile("geocache_"))

    f <- tempfile()
    writeLines("payload", f)
    url <- "https://example.com/GSE1/data.txt.gz"

    expect_null(GEOquery:::.cache_lookup(url))
    GEOquery:::.cache_add(url, f)

    cached <- GEOquery:::.cache_lookup(url)
    expect_false(is.null(cached))
    expect_true(file.exists(cached))
    expect_equal(readLines(cached), "payload")

    clearGEOCache()
    expect_null(GEOquery:::.cache_lookup(url))
})

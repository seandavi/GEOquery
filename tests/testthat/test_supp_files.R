library(testthat)
context('GEO supplementary files')

test_that("GSE Supplemental files downloading works", {
    res = getGEOSuppFiles('GSE1000')
    expect_equivalent(nrow(res), 1)
})

test_that("GSM Supplemental files downloading works", {
    res = getGEOSuppFiles('GSM15789')
    expect_equivalent(nrow(res), 1)
})

test_that("GSM Supplemental files downloading to baseDir works", {
    res = getGEOSuppFiles('GSM15789', baseDir=tempdir())
    expect_equivalent(nrow(res), 1)
})

test_that("GSE supplemental file no-download works", {
    res = getGEOSuppFiles('GSE63137', fetch_files = FALSE)
    expect_equivalent(nrow(res), 12)
    expect_equivalent(ncol(res), 2)
})

test_that("GSE Supplemental file filtering works", {
    res = getGEOSuppFiles('GSE63137', fetch_files = FALSE, filter_regex = 'txt.gz')
    expect_equivalent(ncol(res), 2)
    expect_equivalent(nrow(res), 4)
})

library(testthat)
context('GEO supplementary files')

test_that("GSE Supplemental files downloading works", {
    res = getGEOSuppFiles('GSE1000')
    expect_equivalent(nrow(res),1)
})

test_that("GSM Supplemental files downloading works", {
    res = getGEOSuppFiles('GSM15789')
    expect_equivalent(nrow(res),1)
})

test_that("GSM Supplemental files downloading to baseDir works", {
    res = getGEOSuppFiles('GSM15789', baseDir=tempdir())
    expect_equivalent(nrow(res),1)
})

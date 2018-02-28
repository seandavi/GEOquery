library(GEOquery)
context("GSE")


test_that("empty GSE is handled correctly", {
    gse = getGEO('GSE11413')

    expect_is(gse, 'list')
    expect_is(gse[[1]], 'ExpressionSet')
    expect_equal(nrow(pData(gse[[1]])),12)
    expect_equal(nrow(fData(gse[[1]])),0)
})

test_that("case-mismatched IDs in GSEs handled correctly", {
    gse = getGEO('GSE35683')

    expect_equivalent(nrow(gse[[1]]), 54675)
})

test_that("single-sample GSE handled correctly", {
    gse = getGEO("GSE11595")
    
    expect_is(gse[[1]], "ExpressionSet")
    expect_equivalent(ncol(gse[[1]]),1)
})

test_that("short GSE handled correctly", {
    gse = getGEO("GSE34145")

    expect_equivalent(nrow(gse[[1]]), 15)
})

test_that("generic SOFT format GSE handled correctly", {
    gse = getGEO('GSE1563',GSEMatrix=FALSE)

    expect_equal(62, length(GSMList(gse)))
    expect_equal(1, length(GPLList(gse)))
    expect_equal(12625, nrow(Table(GPLList(gse)[[1]])))
    expect_equal(12625, nrow(Table(GSMList(gse)[[1]])))
    sapply(GSMList(gse),function(x) {
        expect_is(x,'GSM') # all GSMList members should be "GSM"s
        expect_equivalent(12625, nrow(Table(x)))}) # all GSMs should have the 12625 rows
    expect_is(GPLList(gse)[[1]],'GPL')
})

test_that("GSE with more than one value per characteristic handled", {
  gse = getGEO("GSE71989")
  
  expect_equivalent(nrow(gse[[1]]), 54675)
  expect_equivalent(ncol(gse[[1]]), 22)
})

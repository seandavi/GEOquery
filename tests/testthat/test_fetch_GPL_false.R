library(GEOquery)
context('get records without GPL')

test_that("GSE without GPL works", {
    gse = getGEO('GSE2553',getGPL=FALSE)[[1]]
    
    expect_true(validObject(gse))
    expect_equivalent(0,ncol(fData(gse)))
    expect_equivalent(12600,nrow(fData(gse)))
})


test_that("GDS without GPL works", {
    gds = getGEO('GDS10')
    
    eset = GDS2eSet(gds,getGPL=FALSE)
    expect_true(validObject(eset))
    expect_equivalent(0,ncol(fData(eset)))
    expect_equivalent(39114,nrow(fData(eset)))
})

library(testthat)
library(limma)
context('GEO conversions')

test_that("GDS2eSet works", {
    gds = getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))

    eset = GDS2eSet(gds)

    expect_is(eset,'ExpressionSet') # eset should be an ExpressionSet
    expect_equivalent(pubMedIds(experimentData(eset)),'14641932') #basic experimentData check failed
    expect_equivalent(nrow(eset),22645) #'eset has wrong number of rows')
    expect_equivalent(ncol(eset),17) #'eset has wrong number of columns')
})
    
test_that("GDS2MA works", {
    gds = getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))
    
    malist = GDS2MA(gds)

    expect_is(malist,'MAList') #,'malist should be an MAList')
    expect_equivalent(nrow(malist),22645) #'malist has wrong number of rows')
    expect_equivalent(ncol(malist),17) #malist has wrong number of columns')
})

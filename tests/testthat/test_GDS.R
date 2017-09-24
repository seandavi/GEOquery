library(GEOquery)
context('GDS')

test_that("generic GDS parsing works as expected", {
    gds = getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))

    expect_is(gds,'GDS') # gds should be a "GDS"
    expect_is(Meta(gds),'list') # Meta(gds) should be a list
    expect_is(Table(gds),'data.frame') #Table(gds) should be a data.frame
    expect_is(dataTable(gds),'GEODataTable') #dataTable(gds) should be a GEODataTable
    expect_equivalent(nrow(Table(gds)),22645) # gds has wrong number of rows!
    expect_equivalent(ncol(Table(gds)),19) #gds has wrong number of columns!
    expect_length(Meta(gds),23) # gds has wrong number of Metadata entries!
})

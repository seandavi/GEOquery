library(GEOquery)
context('GPL')

test_that("generic GPL parsing works as expected", {
    gpl = getGEO('GPL96')
    
    expect_is(gpl,'GPL') #gpl is not a GPL object!')
    expect_is(gpl,'GEOData') #gpl does not inherit from GEOData!')
    expect_equivalent(nrow(Table(gpl)),22283) # nrow(gpl) does not match!')
    expect_equivalent(ncol(Table(gpl)),16) #ncol(gpl) does not match!')
    expect_is(Meta(gpl), 'list')  #Meta(gpl) should be a list')
    expect_is(Table(gpl), 'data.frame') #Table(gpl) should be a data.frame')
    expect_is(dataTable(gpl), 'GEODataTable') #dataTable(gpl) should be a GEODataTable')
})

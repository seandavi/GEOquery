library(testthat)
context('GSM')

test_that("basic GSM works", {
    gsm = getGEO('GSM11805')

    expect_is(gsm,'GSM')
    expect_is(Meta(gsm),'list')
    expect_is(Table(gsm),'data.frame')
    expect_is(dataTable(gsm),'GEODataTable')
    expect_equivalent(Accession(gsm),'GSM11805')
    expect_equivalent(nrow(Table(gsm)),22283)
    expect_equivalent(ncol(Table(gsm)),3)
    expect_length(Meta(gsm),28)
})

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


test_that("GSE has populated experimentData", {
  gse = getGEO("GSE53986")
  
  ed <- experimentData(gse[[1]])
  expect_equal(pubMedIds(ed), "24739962")
  
  ei <- expinfo(ed)
  expect_equivalent(ei[1], "Jason,A,Hackney")
  expect_equivalent(ei[2], "") #lab
  expect_equivalent(ei[3], "hackney.jason@gene.com")
  expect_equivalent(ei[4], "NRROS negatively regulates ROS in phagocytes during host defense and autoimmunity")
  expect_equivalent(ei[5], "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53986") #url
})

test_that("GSE populates experimentData as much as possible", {
  gse = getGEO("GSE99709")
  
  ed <- experimentData(gse[[1]])
  expect_equal(pubMedIds(ed), "")
  
  ei <- expinfo(ed)
  expect_equivalent(ei[1], "John,,Mariani")
  expect_equivalent(ei[2], "") #lab
  expect_equivalent(ei[3], "john_mariani@urmc.rochester.edu")
  expect_equivalent(ei[4], "RNA-Sequencing of Stat3 silenced oligodendrocyte progenitor cells.")
  expect_equivalent(ei[5], "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99709") #url
  # ----------------------------------------------------------------
  gse = getGEO("GSE27712")
  
  ed <- experimentData(gse[[1]])
  expect_equal(pubMedIds(ed), "22253802")
   
  ei <- expinfo(ed)
  expect_equivalent(ei[1], "Joachim,L,Schultze")
  expect_equivalent(ei[2], "") #lab
  expect_equivalent(ei[3], "j.schultze@uni-bonn.de")
  expect_equivalent(ei[4], "GC424 tumor cells and gastric tumors")
  expect_equivalent(ei[5], "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27712") #url
  expect_equivalent(abstract(ed), "This SuperSeries is composed of the SubSeries listed below.")
})

test_that("Empty files produces an error", {
  suppressWarnings(
    expect_error(getGEO(filename = system.file('extdata/GPLbroken.soft.gz', package="GEOquery")))
  )
})

test_that("GSE/GPL with integer64 columns handled correctly", {
  gse = getGEO("GSE7864")[[1]]
  fdata = fData(gse)
  expect_s3_class(fdata$ID, "integer64")
  expect_is(rownames(fdata), "character")
})

library(GEOquery)

#############################
#
# GPL tests
#
#############################
context('GPL tests')

test_that('Short GPL files are treated correctly',{
  gpl = getGEO('GPL15505')
  expect_that(inherits(gpl,'GPL'),   is_true())
  expect_that(nrow(Table(gpl)),      equals(52))
})

test_that('getGEO on a GPL works correctly',{
  gpl = getGEO('GPL96')

  expect_that(gpl,               is_a('GPL'))
  expect_that(gpl,               is_a('GEOData'))
  expect_that(Table(gpl),        is_a('data.frame'))
  expect_that(nrow(Table(gpl)),  equals(22283))
  expect_that(ncol(Table(gpl)),  equals(16))
  expect_that(Meta(gpl),         is_a('list'))
  expect_that(dataTable(gpl),    is_a('GEODataTable'))
})

#############################
#
# GSEMatrix tests
#
#############################
context('GSEMatrix tests')

test_that('Single Sample GSE works',{
  gse = getGEO('GSE11595')

  expect_is(gse[[1]],"ExpressionSet")
  expect_equivalent(ncol(gse[[1]]),1)
})

test_that('Empty GSE handled correctly',{
  gse <- getGEO('GSE11413')

  expect_that(gse,                   is_a('list'))
  expect_that(gse[[1]],              is_a('ExpressionSet'))
  expect_that(nrow(pData(gse[[1]])), equals(12))
  expect_that(nrow(fData(gse[[1]])), equals(0))
})

test_that('GSE with quoted IDREF works correctly',{
  gse = getGEO('GSE19711')

  expect_that(length(gse),    equals(3))
})

test_that('GSE with case-mismatched IDs works',{
  gse = getGEO('GSE35683')

  expect_that(nrow(gse[[1]]),       is_equivalent_to(54675))
})

#############################
#
# GSE SOFT format tests
#
#############################

context('GSE SOFT format tests')

test_that('GSE SOFT format works',{
  gse = getGEO('GSE1563',GSEMatrix=FALSE)

  expect_that(length(GSMList(gse)),      equals(62))
  expect_that(length(GPLList(gse)),      equals(1))
  expect_that(nrow(Table(GPLList(gse)[[1]])), equals(12625))
  expect_that(nrow(Table(GSMList(gse)[[1]])), equals(12625))
  sapply(GSMList(gse),function(x) {
    expect_that(x,                       is_a('GSM'))})
  expect_that(GPLList(gse)[[1]],         is_a('GPL'))
})

#############################
#
# GDS tests
#
#############################

context('GDS tests')

test_that('Basic GDS functionality',{
  gds = getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))

  expect_is(gds,'GDS')
  expect_is(Meta(gds),'list')
  expect_is(Table(gds),'data.frame')
  expect_is(dataTable(gds),'GEODataTable')
  expect_equal(nrow(Table(gds)),22645)
  expect_equal(ncol(Table(gds)),19)
  expect_equal(length(Meta(gds)),23)
})

test_that('GDS conversions',{
  gds = getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))
  eset = GDS2eSet(gds)
  
  expect_is(eset,'ExpressionSet')
  expect_equal(pubMedIds(experimentData(eset)),'14641932')
  expect_equivalent(nrow(eset),22645)
  expect_equivalent(ncol(eset),17)

  malist = GDS2MA(gds)

  expect_is(malist,'MAList')
  expect_equivalent(nrow(malist),22645)
  expect_equivalent(ncol(malist),17)
})
  


#############################
#
# GSM tests
#
#############################

context('GSM tests')

test_that('GSM functionality',{
  gsm = getGEO(filename=system.file("extdata/GSM11805.txt.gz",package="GEOquery"))

  expect_is(gsm,'GSM')
  expect_is(Meta(gsm),'list')
  expect_is(Table(gsm),'data.frame')
  expect_is(dataTable(gsm),'GEODataTable')
  expect_equal(Accession(gsm),'GSM11805')
  expect_equal(nrow(Table(gsm)),22283)
  expect_equal(ncol(Table(gsm)),3)
  expect_equal(length(Meta(gsm)),28)
})

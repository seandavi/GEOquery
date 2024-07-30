library(testthat)
context("geo_rnaseq")

test_that("getGSEDownloadPageURLs returns character vector", {
  gse <- "GSE83322"
  result <- getGSEDownloadPageURLs(gse)

  expect_type(result, "character")
  expect_true(length(result) > 0)
})

test_that("getRNAQuantRawCountsURL returns character vector", {
  gse <- "GSE83322"
  links <- getGSEDownloadPageURLs(gse)
  result <- getRNAQuantRawCountsURL(links)

  expect_type(result, "character")
  expect_length(result, 1)
})

test_that("getRNAQuantAnnotationURL returns character vector", {
  gse <- "GSE83322"
  links <- getGSEDownloadPageURLs(gse)
  result <- getRNAQuantAnnotationURL(links)

  expect_type(result, "character")
  expect_length(result, 1)
})

links <- getGSEDownloadPageURLs("GSE83322")

test_that("readRNAQuantRawCounts returns matrix", {
  link <- getRNAQuantRawCountsURL(links)
  result <- readRNAQuantRawCounts(link)

  expect_true(is.matrix(result))
  expect_true(nrow(result) > 0)
  expect_true(ncol(result) > 0)
})

test_that("readRNAQuantAnnotation returns data.frame", {
  link <- getRNAQuantAnnotationURL(links)
  result <- readRNAQuantAnnotation(link)

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true(ncol(result) > 0)
})

test_that("getRNASeqQuantResults returns list", {
  gse <- "GSE83322"
  result <- getRNASeqQuantResults(gse)

  expect_type(result, "list")
  expect_length(result, 2)
  expect_true(is.matrix(result$quants))
  expect_true(is.data.frame(result$annotation))
})

test_that("getRNASeqData returns SummarizedExperiment", {
  gse <- "GSE83322"
  result <- getRNASeqData(gse)

  expect_s4_class(result, "SummarizedExperiment")
  expect_true(nrow(result) > 0)
  expect_true(ncol(result) > 0)
})

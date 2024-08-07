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

test_that("hasRNASeqQuantifications returns logical", {
  gse <- "GSE83322"
  result <- hasRNASeqQuantifications(gse)

  expect_type(result, "logical")
})

test_that("hasRNASeqQuantifications returns TRUE for GSE83322", {
  gse <- "GSE83322"
  result <- hasRNASeqQuantifications(gse)

  expect_true(result)
})

test_that("hasRNASeqQuantifications returns FALSE for GSE2553", {
  gse <- "GSE2553"
  result <- hasRNASeqQuantifications(gse)

  expect_false(result)
})

test_that("getRNASeqQuantGenomeInfo returns correct data", {
  genome_info <- getRNASeqQuantGenomeInfo("GSE83322")

  expect_length(genome_info, 3)
  expect_true(genome_info["species"] == "Human")
  expect_true(stringr::str_starts(genome_info["genome_build"], "GR"))
})

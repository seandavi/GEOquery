library(testthat)
context("geo_rnaseq")

test_that("test_returns_character_vector", {
  gse <- "GSE83322"
  result <- getGSEDownloadPageURLs(gse)

  expect_type(result, "character")
  expect_true(length(result) > 0)
})

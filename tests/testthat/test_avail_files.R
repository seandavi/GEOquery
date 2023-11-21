library(GEOquery)
context("gse_download_types")

test_that("GSE download file manifest works", {
  res <- available_gse_files("GSE83322")
  descriptions <- dplyr::pull(res, description)
  expect_gte(nrow(res), 7)
})

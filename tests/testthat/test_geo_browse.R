context("browsing geo with accession number")

test_that("returns_correct_url_for_valid_geo_accession", {
  geo <- "GSE262484"
  expected_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262484"
  result <- urlForAccession(geo)
  expect_equal(result, expected_url)
})

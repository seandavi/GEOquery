# Offline unit tests for the JSON esummary backend of searchGEO() (R/searchGEO.R,
# #223). The network fetch (.fetch_gds_esummary_json) is integration-only; here
# we test the pure parser against a saved gds esummary (v2.0) JSON fixture.

test_that(".parse_gds_esummary_json() parses records into the tidy table", {
  json <- paste(
    readLines(test_path("fixtures", "gds_esummary.json")),
    collapse = "\n"
  )
  res <- GEOquery:::.parse_gds_esummary_json(json)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2L)
  expect_true(all(
    c("Accession", "Title", "Summary", "Organism", "Type", "GPL",
      "n_samples", "PDAT", "suppFile", "FTPLink", "entryType", "ID") %in%
      names(res)
  ))

  expect_identical(res$Accession, c("GSE12345", "GSE67890"))
  expect_identical(res$Title, c("Study A title", "Study B title"))
  expect_identical(res$Organism, c("Homo sapiens", "Mus musculus"))
  # n_samples is coerced to integer.
  expect_identical(res$n_samples, c(20L, 6L))
  expect_identical(res$ID, c("200012345", "200067890"))
  expect_identical(res$GPL, c("570", "16791"))
})

test_that(".parse_gds_esummary_json() returns an empty typed frame for no hits", {
  empty <- '{"header":{},"result":{"uids":[]}}'
  res <- GEOquery:::.parse_gds_esummary_json(empty)
  expect_equal(nrow(res), 0L)
  expect_true(all(c("Accession", "Title", "ID") %in% names(res)))
})

test_that(".gds_summary_row() fills missing fields with NA", {
  row <- GEOquery:::.gds_summary_row(list(accession = "GSE1"))
  expect_equal(nrow(row), 1L)
  expect_identical(row$Accession, "GSE1")
  expect_true(is.na(row$Title))
  expect_true(is.na(row$n_samples))
})

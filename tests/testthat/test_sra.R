# Offline unit tests for the SRA resolver's pure parser (R/sra.R, #219).
# The network wrapper geoToSRA() is exercised only in the integration suite;
# here we test .parse_sra_runinfo() against a saved runinfo CSV fixture.

test_that(".parse_sra_runinfo() maps runinfo columns to the tidy SRA table", {
  runinfo <- paste(
    readLines(test_path("fixtures", "sra_runinfo.csv")),
    collapse = "\n"
  )
  res <- GEOquery:::.parse_sra_runinfo(runinfo, "GSE300001")

  expect_s3_class(res, "data.frame")
  expect_named(res, c("geo_accession", "srp", "srx", "srr", "srs"))
  expect_equal(nrow(res), 3L)

  expect_true(all(res$geo_accession == "GSE300001"))
  # One study across all runs.
  expect_identical(unique(res$srp), "SRP300001")
  # Two experiments, three runs, two samples.
  expect_identical(res$srr, c("SRR12345678", "SRR12345679", "SRR12345680"))
  expect_identical(unique(res$srx), c("SRX9999901", "SRX9999902"))
  expect_identical(unique(res$srs), c("SRS8000001", "SRS8000002"))
})

test_that(".parse_sra_runinfo() returns an empty table for a header-only body", {
  header <- readLines(test_path("fixtures", "sra_runinfo.csv"))[1L]
  res <- GEOquery:::.parse_sra_runinfo(header, "GSE300001")
  expect_equal(nrow(res), 0L)
  expect_named(res, c("geo_accession", "srp", "srx", "srr", "srs"))
})

test_that(".parse_sra_runinfo() returns an empty table for unrelated CSV", {
  res <- GEOquery:::.parse_sra_runinfo("a,b\n1,2\n", "GSE300001")
  expect_equal(nrow(res), 0L)
})

test_that("geoToSRA() rejects non-GEO input before any network call", {
  expect_error(GEOquery::geoToSRA("not-an-accession"), class = "geoquery_bad_accession")
  expect_error(GEOquery::geoToSRA(c("GSE1", "GSE2")), class = "geoquery_bad_accession")
  expect_error(GEOquery::geoToSRA("GPL570"), class = "geoquery_bad_accession")
})

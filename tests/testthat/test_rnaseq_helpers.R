# Offline unit tests for the pure RNA-seq URL/link helpers in R/rnaseq.R.
# These construct all inputs directly and never touch the network.

annot_url <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/download/",
  "?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
)

test_that("extractFilenameFromDownloadURL() returns the file query value", {
  expect_equal(
    GEOquery:::extractFilenameFromDownloadURL(annot_url),
    "Human.GRCh38.p13.annot.tsv.gz"
  )
})

test_that("extractFilenameFromDownloadURL() returns NULL without a file param", {
  no_file_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164073"
  expect_null(GEOquery:::extractFilenameFromDownloadURL(no_file_url))
})

test_that("extractFilenameFromDownloadURL() returns NULL for character(0)", {
  # Regression test: previously errored because httr2::url_parse() requires
  # a single string.
  expect_null(GEOquery:::extractFilenameFromDownloadURL(character(0)))
})

test_that(".extractGeoDownloadLinks() parses and normalizes hrefs offline", {
  html <- paste0(
    "<html><body>",
    "<a href='/geo/download/?acc=GSE1&file=raw_counts.tsv.gz'>raw</a>",
    "<a href='ftp://ftp.ncbi.nlm.nih.gov/geo/GSE1_annot.tsv.gz'>annot</a>",
    "<a>no-href</a>",
    "</body></html>"
  )
  links <- GEOquery:::.extractGeoDownloadLinks(html)

  expect_s3_class(links, "geoDownloadLinks")
  # leading /geo/ and ftp:// are rewritten to absolute https URLs
  expect_true(any(grepl(
    "^https://www.ncbi.nlm.nih.gov/geo/download/\\?acc=GSE1&file=raw_counts",
    links
  )))
  expect_true(any(grepl(
    "^https://ftp.ncbi.nlm.nih.gov/geo/GSE1_annot.tsv.gz$", links
  )))
  # the raw-counts / annotation selectors still work on the result
  expect_match(GEOquery:::getRNAQuantRawCountsURL(links), "raw_counts")
  expect_match(GEOquery:::getRNAQuantAnnotationURL(links), "annot.tsv.gz")
})

test_that("urlExtractRNASeqQuantGenomeInfo() parses genome build and species", {
  info <- GEOquery:::urlExtractRNASeqQuantGenomeInfo(annot_url)
  expect_equal(info[["genome_build"]], "GRCh38.p13")
  expect_equal(info[["species"]], "Human")
  expect_equal(info[["fname"]], "Human.GRCh38.p13.annot.tsv.gz")
})

test_that("urlExtractRNASeqQuantGenomeInfo() returns NULL for missing input", {
  expect_null(GEOquery:::urlExtractRNASeqQuantGenomeInfo(character(0)))
  no_file_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164073"
  expect_null(GEOquery:::urlExtractRNASeqQuantGenomeInfo(no_file_url))
})

test_that("getRNAQuantRawCountsURL()/getRNAQuantAnnotationURL() pick correct URL", {
  raw_url <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/download/",
    "?type=rnaseq_counts&acc=GSE164073&format=file&file=GSE164073_raw_counts",
    "_GRCh38.p13_NCBI.tsv.gz"
  )
  other_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164073"
  links <- structure(
    c(raw_url, annot_url, other_url),
    class = c("geoDownloadLinks", "character")
  )

  expect_equal(GEOquery:::getRNAQuantRawCountsURL(links), raw_url)
  expect_equal(GEOquery:::getRNAQuantAnnotationURL(links), annot_url)
})

test_that("raw counts/annotation helpers require a geoDownloadLinks object", {
  plain <- c("https://example.com/raw_counts", "https://example.com/annot.tsv.gz")
  expect_error(
    GEOquery:::getRNAQuantRawCountsURL(plain),
    "Input must be a geoDownloadLinks object"
  )
  expect_error(
    GEOquery:::getRNAQuantAnnotationURL(plain),
    "Input must be a geoDownloadLinks object"
  )
})

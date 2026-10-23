test_that("readRNAQuantAnnotation errors clearly on a non-annotation file (#235)", {
  # NCBI can serve an HTML reCAPTCHA page instead of the TSV
  html <- tempfile(fileext = ".tsv")
  writeLines(c("<!doctype html><html>", "<body>captcha</body>"), html)
  expect_error(GEOquery:::readRNAQuantAnnotation(html), "no 'GeneID' column")
})

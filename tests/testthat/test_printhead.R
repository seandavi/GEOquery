# Offline test for the internal printHead() helper used by the GEODataTable
# show method. Regression guard: a corrupted definition once made printHead()
# return a function instead of printing the leading rows/elements.

test_that("printHead() prints the leading 5 elements of a long vector", {
  out <- paste(capture.output(GEOquery:::printHead(1:100)), collapse = " ")
  expect_match(out, "1 2 3 4 5")
  expect_match(out, "95 more elements")
})

test_that("printHead() prints a short vector in full", {
  out <- paste(capture.output(GEOquery:::printHead(1:3)), collapse = " ")
  expect_match(out, "1 2 3")
  expect_no_match(out, "more elements")
})

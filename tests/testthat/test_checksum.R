# Offline unit tests for the checksum helpers in R/checksum.R (#222).
# Fully offline: a fixture file is written to a temp path and hashed.

# Known MD5 of the exact bytes written below (computed with tools::md5sum).
.fixture_bytes <- "GEOquery checksum fixture\n"
.fixture_md5 <- "cfe23436c7c6cd4c5d69188b1fdd23d6"

write_fixture <- function() {
  p <- tempfile(fileext = ".txt")
  # writeBin to control the exact bytes (no platform newline translation).
  writeBin(charToRaw(.fixture_bytes), p)
  p
}

test_that(".file_md5() returns the plain lowercase hex digest", {
  p <- write_fixture()
  on.exit(unlink(p))
  expect_identical(GEOquery:::.file_md5(p), .fixture_md5)
})

test_that(".file_md5() errors on a missing file", {
  expect_error(
    GEOquery:::.file_md5(tempfile()),
    class = "geoquery_download_error"
  )
})

test_that(".verify_md5() returns TRUE on a match, without warning", {
  p <- write_fixture()
  on.exit(unlink(p))
  expect_silent(res <- GEOquery:::.verify_md5(p, .fixture_md5))
  expect_true(res)
  # Case-insensitive comparison.
  expect_true(GEOquery:::.verify_md5(p, toupper(.fixture_md5)))
})

test_that(".verify_md5() warns (structured) and returns FALSE on mismatch", {
  p <- write_fixture()
  on.exit(unlink(p))
  expect_warning(
    res <- GEOquery:::.verify_md5(p, "00000000000000000000000000000000"),
    class = "geoquery_checksum_mismatch"
  )
  expect_false(res)
})

test_that(".verify_md5() is a no-op when no expected value is supplied", {
  p <- write_fixture()
  on.exit(unlink(p))
  expect_true(GEOquery:::.verify_md5(p, NULL))
  expect_true(GEOquery:::.verify_md5(p, NA_character_))
  expect_true(GEOquery:::.verify_md5(p, character(0)))
})

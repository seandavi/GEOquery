# Checksum verification for downloaded / cached files (#222).
#
# A truncated or corrupted download, or a stale cache entry, is otherwise
# trusted silently. These small internal helpers let the download path
# optionally verify a file against an expected MD5 (e.g. one supplied by the
# ENA raw-file resolver (#220) or recorded in a fetch receipt (#221)).
# Verification is advisory: a mismatch warns via a structured condition rather
# than erroring, so callers decide whether to retry or discard.

# Compute a file's MD5 as a plain lowercase hex string (tools::md5sum returns a
# vector named by path; drop the name).
.file_md5 <- function(path) {
  if (!file.exists(path)) {
    .abort_download(sprintf("Cannot checksum missing file '%s'.", path))
  }
  unname(tools::md5sum(path))
}

# Verify a file's MD5 against `expected`. Returns TRUE on match (or when there
# is nothing to check, i.e. `expected` is NULL/NA). On mismatch emits a
# `geoquery_checksum_mismatch` warning and returns FALSE. Comparison is
# case-insensitive (hex digests may arrive in either case).
.verify_md5 <- function(path, expected) {
  if (is.null(expected) || length(expected) == 0L || is.na(expected)) {
    return(TRUE)
  }
  actual <- .file_md5(path)
  if (!identical(tolower(actual), tolower(as.character(expected)))) {
    .warn_checksum_mismatch(path = path, expected = expected, actual = actual)
    return(FALSE)
  }
  TRUE
}

# Offline unit tests for the pure helpers in R/searchGEO.R.
#
# These functions are all internal (non-exported) and are reached with the
# `GEOquery:::` triple-colon operator. Every input is constructed directly in
# the test -- nothing here touches the network (no rentrez / Entrez calls).

# ---------------------------------------------------------------------------
# str_split(): split each string ONCE on the first match of a perl regex,
# returning a list of length-2 (before, after) character vectors.
# ---------------------------------------------------------------------------

test_that("str_split() splits a single string once into (before, after)", {
  res <- GEOquery:::str_split("age: 30", "(\\s*+):(\\s*+)")
  expect_type(res, "list")
  expect_length(res, 1L)
  expect_identical(res[[1L]], c("age", "30"))
})

test_that("str_split() is vectorized over the input strings", {
  res <- GEOquery:::str_split(c("age: 30", "sex: F"), "(\\s*+):(\\s*+)")
  expect_type(res, "list")
  expect_length(res, 2L)
  expect_identical(res[[1L]], c("age", "30"))
  expect_identical(res[[2L]], c("sex", "F"))
})

test_that("str_split() splits only on the FIRST match", {
  # Two colons: only the first should be used as the split point, so the
  # trailing colon stays inside the second element.
  res <- GEOquery:::str_split("a:b:c", ":")
  expect_identical(res[[1L]], c("a", "b:c"))
})

test_that("str_split() returns the whole string when there is no match", {
  res <- GEOquery:::str_split("nomatch", "(\\s*+):(\\s*+)")
  expect_length(res, 1L)
  expect_identical(res[[1L]], "nomatch")
})

test_that("str_split() honors ignore.case", {
  res <- GEOquery:::str_split("KEY: v", "key", ignore.case = TRUE)
  expect_identical(res[[1L]], c("", ": v"))

  # Without ignore.case the lowercase pattern must not match the uppercase key.
  res_cs <- GEOquery:::str_split("KEY: v", "key")
  expect_identical(res_cs[[1L]], "KEY: v")
})

# ---------------------------------------------------------------------------
# parse_name_value_pairs(): turn a list of "key:value" character vectors into
# a list keyed by unique key names, each element the vector of values across
# the records (filling missing keys across records).
# ---------------------------------------------------------------------------

test_that("parse_name_value_pairs() keys by name and collects values per record", {
  res <- GEOquery:::parse_name_value_pairs(
    list(c("age:30", "sex:F"), c("age:40", "sex:M"))
  )
  expect_type(res, "list")
  expect_named(res, c("age", "sex"))
  # Numeric-looking values are coerced to numeric via data.table::fread.
  expect_equal(res$age, c(30, 40))
  expect_identical(res$sex, c("F", "M"))
})

test_that("parse_name_value_pairs() strips whitespace around the separator", {
  res <- GEOquery:::parse_name_value_pairs(list(c("age: 30", "sex: F")))
  expect_named(res, c("age", "sex"))
  expect_equal(res$age, 30)
  expect_identical(res$sex, "F")
})

test_that("parse_name_value_pairs() fills missing keys across records with NA", {
  res <- GEOquery:::parse_name_value_pairs(
    list(c("age:30", "sex:F"), c("age:40", "weight:70"))
  )
  expect_named(res, c("age", "sex", "weight"))
  expect_equal(res$age, c(30, 40))
  # Record 2 has no "sex"; record 1 has no "weight" -> filled with NA.
  expect_identical(res$sex, c("F", NA))
  expect_equal(res$weight, c(NA, 70))
})

test_that("parse_name_value_pairs() handles an empty record", {
  res <- GEOquery:::parse_name_value_pairs(list(character(0)))
  expect_type(res, "list")
  expect_length(res, 0L)
})

# ---------------------------------------------------------------------------
# preprocess_records(): normalize a single Entrez esummary-style text block
# into a list (one element per record) of relabeled, newline-split lines.
# ---------------------------------------------------------------------------

test_that("preprocess_records() relabels an esummary record into keyed lines", {
  rec <- paste(
    "1. My great study title",
    "(Submitter supplied) This is the summary of the study.",
    "Organism:\tHomo sapiens",
    "Type:\t\tExpression profiling by array",
    "Platform: GPL570 20 Samples",
    "Series\t\tAccession: GSE12345\tID: 200012345",
    sep = "\n"
  )
  out <- GEOquery:::preprocess_records(rec)

  expect_type(out, "list")
  expect_length(out, 1L)

  lines <- out[[1L]]
  # Leading "N." numbering becomes a "Title:" prefix.
  expect_true("Title: My great study title" %in% lines)
  # "(Submitter supplied)" becomes a "Summary:" line.
  expect_true("Summary: This is the summary of the study." %in% lines)
  # "Platform:" is pluralized to "Platforms:".
  expect_true("Platforms: GPL570" %in% lines)
  # The tab-prefixed "ID:" is pulled onto its own "ID:" line.
  expect_true("ID: 200012345" %in% lines)
  # Tab-collapsed key/value lines survive as clean "key: value".
  expect_true("Organism: Homo sapiens" %in% lines)
})

test_that("preprocess_records() feeds cleanly into parse_name_value_pairs()", {
  rec <- paste(
    "1. Title A",
    "Organism:\tHomo sapiens",
    "Series\t\tAccession: GSE1\tID: 200000001",
    sep = "\n"
  )
  parsed <- GEOquery:::parse_name_value_pairs(
    GEOquery:::preprocess_records(rec)
  )
  expect_true("Title" %in% names(parsed))
  expect_true("Organism" %in% names(parsed))
  expect_true("ID" %in% names(parsed))
  expect_identical(parsed$Title, "Title A")
  expect_identical(parsed$Organism, "Homo sapiens")
  expect_equal(parsed$ID, 200000001)
})

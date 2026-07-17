test_that("txtGrab returns the matching substring", {
    expect_equal(GEOquery:::txtGrab("[0-9]+", "abc123def"), "123")
    expect_equal(GEOquery:::txtGrab("^abc", "abcdef"), "abc")
})

test_that("txtGrab returns empty string when there is no match", {
    expect_equal(GEOquery:::txtGrab("[0-9]+", "abcdef"), "")
})

test_that("txtGrab coerces non-character input via as.character", {
    # numeric input is coerced to its character representation
    expect_equal(GEOquery:::txtGrab("[0-9]+", 12345), "12345")
    expect_equal(GEOquery:::txtGrab("[0-9]", 12345), "1")
    # factor input is coerced to its character label (not integer codes)
    f <- factor("group7")
    expect_equal(GEOquery:::txtGrab("[0-9]+", f), "7")
})

test_that("txtGrab is vectorized over a character vector", {
    x <- c("abc123", "def456", "nomatch")
    expect_equal(GEOquery:::txtGrab("[0-9]+", x), c("123", "456", ""))
})

test_that("filegrep finds matching line numbers and text", {
    tmp <- tempfile()
    writeLines(
        c(
            "line one no match",
            "^SAMPLE = GSM1",
            "some data",
            "^SAMPLE = GSM2",
            "more data"
        ),
        tmp
    )
    con <- file(tmp, "r")
    on.exit({
        close(con)
        unlink(tmp)
    })

    res <- GEOquery:::filegrep(con, "^\\^SAMPLE")

    expect_s3_class(res, "data.frame")
    expect_named(res, c("foundLines", "foundTypes"))
    expect_equal(res$foundLines, c(2, 4))
    expect_equal(res$foundTypes, c("^SAMPLE = GSM1", "^SAMPLE = GSM2"))
})

test_that("filegrep line numbers are correct across chunk boundaries", {
    tmp <- tempfile()
    writeLines(
        c(
            "nomatch 1",
            "MATCH a",
            "nomatch 2",
            "nomatch 3",
            "MATCH b",
            "nomatch 4",
            "MATCH c"
        ),
        tmp
    )
    con <- file(tmp, "r")
    on.exit({
        close(con)
        unlink(tmp)
    })

    # chunksize = 2 forces several readLines() chunks; global line numbers
    # must remain correct across the chunk boundaries.
    res <- GEOquery:::filegrep(con, "^MATCH", chunksize = 2)

    expect_equal(res$foundLines, c(2, 5, 7))
    expect_equal(res$foundTypes, c("MATCH a", "MATCH b", "MATCH c"))
})

test_that("filegrep returns NULL when there are no matches", {
    tmp <- tempfile()
    writeLines(c("alpha", "beta", "gamma"), tmp)
    con <- file(tmp, "r")
    on.exit({
        close(con)
        unlink(tmp)
    })

    res <- GEOquery:::filegrep(con, "no-such-pattern")
    expect_null(res)
})

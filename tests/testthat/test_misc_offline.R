# Offline unit tests for pure helper functions -- no network access.

# ---------------------------------------------------------------------------
# fileOpen(): opens a connection, handling plain and gzipped files.
# ---------------------------------------------------------------------------
test_that("fileOpen() reads a plain file", {
    content <- c("line one", "line two", "line three")
    tf <- tempfile(fileext = ".txt")
    writeLines(content, tf)
    on.exit(unlink(tf), add = TRUE)

    con <- GEOquery:::fileOpen(tf)
    expect_true(inherits(con, "connection"))
    on.exit(close(con), add = TRUE)
    expect_equal(readLines(con), content)
})

test_that("fileOpen() reads a gzipped file via gzfile", {
    content <- c("gz line one", "gz line two")
    tf <- tempfile(fileext = ".gz")
    gz <- gzfile(tf, "w")
    writeLines(content, gz)
    close(gz)
    on.exit(unlink(tf), add = TRUE)

    con <- GEOquery:::fileOpen(tf)
    expect_true(inherits(con, "connection"))
    on.exit(close(con), add = TRUE)
    expect_equal(readLines(con), content)
})

test_that("fileOpen() errors on a nonexistent file", {
    missing <- tempfile(fileext = ".txt")
    expect_error(GEOquery:::fileOpen(missing), "does not appear to exist")
})

# ---------------------------------------------------------------------------
# getGEOSuppFileURL(): pure FTP suppl/ directory URL builder.
# ---------------------------------------------------------------------------
test_that("getGEOSuppFileURL() builds correct FTP URLs", {
    expect_equal(
        getGEOSuppFileURL("GSM12345"),
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM12nnn/GSM12345/suppl/"
    )
    expect_equal(
        getGEOSuppFileURL("GSE1000"),
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1000/suppl/"
    )
    expect_equal(
        getGEOSuppFileURL("GPL570"),
        "https://ftp.ncbi.nlm.nih.gov/geo/platform/GPLnnn/GPL570/suppl/"
    )
})

# ---------------------------------------------------------------------------
# urlForAccession() / browseGEOAccession() / browseWebsiteRNASeqSearch():
# never open a real browser -- stub browseURL in the GEOquery namespace.
# ---------------------------------------------------------------------------
test_that("urlForAccession() builds the GEO accession page URL", {
    expect_equal(
        urlForAccession("GSE262484"),
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262484"
    )
})

test_that("browseGEOAccession() passes urlForAccession()'s URL to browseURL", {
    recorded <- NULL
    testthat::local_mocked_bindings(
        browseURL = function(url, ...) {
            recorded <<- url
            invisible(url)
        },
        .package = "utils"
    )
    browseGEOAccession("GSE262484")
    expect_equal(recorded, urlForAccession("GSE262484"))
})

test_that("browseWebsiteRNASeqSearch() passes an rnaseq URL to browseURL", {
    recorded <- NULL
    testthat::local_mocked_bindings(
        browseURL = function(url, ...) {
            recorded <<- url
            invisible(url)
        },
        .package = "utils"
    )
    browseWebsiteRNASeqSearch()
    expect_true(grepl("rnaseq", recorded, ignore.case = TRUE))
})

# ---------------------------------------------------------------------------
# GDS2MA(): offline path with GPL = NULL, getGPL = FALSE (no network fetch).
# ---------------------------------------------------------------------------
make_fake_gds <- function() {
    tbl <- data.frame(
        ID_REF = c("probe_1", "probe_2", "probe_3"),
        IDENTIFIER = c("geneA", "geneB", "geneC"),
        GSM1 = c(1.0, 2.0, 3.0),
        GSM2 = c(4.0, 5.0, 6.0),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    cols <- data.frame(sample = c("GSM1", "GSM2"), stringsAsFactors = FALSE)
    new("GDS",
        header = list(title = "synthetic GDS", description = "fixture",
            platform = "GPL000001"),
        dataTable = new("GEODataTable", columns = cols, table = tbl))
}

test_that("GDS2MA() builds an MAList offline (GPL = NULL, getGPL = FALSE)", {
    skip_if_not_installed("limma")
    gds <- make_fake_gds()
    expect_warning(
        ma <- GDS2MA(gds, GPL = NULL, getGPL = FALSE),
        "deprecated"
    )

    expect_s4_class(ma, "MAList")
    # M holds only the GSM columns as a numeric matrix.
    expect_true(is.matrix(ma$M))
    expect_equal(dim(ma$M), c(3L, 2L))
    expect_equal(colnames(ma$M), c("GSM1", "GSM2"))
    expect_true(is.numeric(ma$M))
    expect_equal(as.numeric(ma$M[, "GSM1"]), c(1, 2, 3))
    expect_equal(as.numeric(ma$M[, "GSM2"]), c(4, 5, 6))
    # No A matrix and no genes on the offline path.
    expect_null(ma$A)
    expect_null(ma$genes)
})

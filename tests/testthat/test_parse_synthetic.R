# Offline, synthetic-fixture tests for the series-matrix parser.
#
# These exercise the parsing logic without any network access or stored data
# files: the "fixture" is a few KB of crafted text written to a tempfile (see
# issue #169 for the testing strategy). Keep fixtures minimal.

# Build a minimal GEO series-matrix file and return its path. The file has
# `n_features` probes x 2 samples, two characteristics_ch1 key:value rows, and
# the table delimited by the !series_matrix_table_begin/_end markers.
make_fake_series_matrix <- function(n_features = 2) {
    probes <- sprintf("probe_%d", seq_len(n_features))
    header <- c(
        "!Series_title\tGEOquery synthetic test series",
        "!Series_geo_accession\tGSE000001",
        "!Series_summary\tsynthetic fixture for offline parser tests",
        "!Sample_title\tsampleA\tsampleB",
        "!Sample_geo_accession\tGSM000001\tGSM000002",
        "!Sample_characteristics_ch1\ttissue: liver\ttissue: brain",
        "!Sample_characteristics_ch1\tagent: none\tagent: drug",
        "!Sample_platform_id\tGPL000001\tGPL000001",
        "!series_matrix_table_begin",
        "\"ID_REF\"\t\"GSM000001\"\t\"GSM000002\""
    )
    rows <- vapply(seq_len(n_features), function(i) {
        sprintf("\"%s\"\t%f\t%f", probes[i], i + 0.5, i + 10.5)
    }, character(1))
    txt <- c(header, rows, "!series_matrix_table_end")
    f <- tempfile(fileext = ".txt")
    writeLines(txt, f)
    f
}

test_that("parseGSEMatrix parses a synthetic series matrix offline (getGPL = FALSE)", {
    f <- make_fake_series_matrix(n_features = 3)
    res <- GEOquery:::parseGSEMatrix(f, getGPL = FALSE)
    eset <- res$eset

    expect_s4_class(eset, "ExpressionSet")
    expect_equal(nrow(Biobase::exprs(eset)), 3L)
    expect_equal(ncol(Biobase::exprs(eset)), 2L)
    expect_equal(Biobase::sampleNames(eset), c("GSM000001", "GSM000002"))
    # feature names come from featureData (probes); assert value by position
    # because exprs() rownames are dropped on the getGPL = FALSE path (a known
    # inconsistency in parseGSEMatrix -- see #173).
    expect_equal(Biobase::featureNames(eset), c("probe_1", "probe_2", "probe_3"))
    expect_equal(unname(Biobase::exprs(eset)[1, 1]), 1.5)
})

test_that("parseCharacteristics = TRUE unpacks characteristics_ch1 into pData columns", {
    f <- make_fake_series_matrix()
    eset <- GEOquery:::parseGSEMatrix(f, getGPL = FALSE, parseCharacteristics = TRUE)$eset
    pd <- Biobase::pData(eset)

    expect_true("tissue:ch1" %in% colnames(pd))
    expect_true("agent:ch1" %in% colnames(pd))
    expect_equal(trimws(as.character(pd["GSM000001", "tissue:ch1"])), "liver")
    expect_equal(trimws(as.character(pd["GSM000002", "agent:ch1"])), "drug")
})

test_that("parseCharacteristics = FALSE leaves characteristics unparsed (regression for #60)", {
    f <- make_fake_series_matrix()
    eset <- GEOquery:::parseGSEMatrix(f, getGPL = FALSE, parseCharacteristics = FALSE)$eset
    cols <- colnames(Biobase::pData(eset))

    # the split key:value columns must NOT be created ...
    expect_false("tissue:ch1" %in% cols)
    expect_false("agent:ch1" %in% cols)
    # ... and the raw characteristics column is retained
    expect_true(any(grepl("^characteristics_ch1", cols)))
})

test_that("parseGEO threads parseCharacteristics through to parseGSEMatrix (#60)", {
    # The #60 bug was that parseCharacteristics was accepted at the top level
    # but dropped before reaching parseGSEMatrix. Exercise the local-file
    # forwarding path (parseGEO -> parseGSEMatrix) offline via getGPL = FALSE.
    f <- make_fake_series_matrix()

    eset_off <- GEOquery:::parseGEO(f, GSElimits = NULL, getGPL = FALSE,
        parseCharacteristics = FALSE)
    expect_false("tissue:ch1" %in% colnames(Biobase::pData(eset_off)))

    eset_on <- GEOquery:::parseGEO(f, GSElimits = NULL, getGPL = FALSE,
        parseCharacteristics = TRUE)
    expect_true("tissue:ch1" %in% colnames(Biobase::pData(eset_on)))
})

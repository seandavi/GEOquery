# Offline regression tests for the issue-triage fixes (#98, #14, #80, #148).
# All fixtures are crafted in-test; no network access (see #169 strategy).

# --- #98: duplicate feature IDs in a series matrix -------------------------

# Build a series matrix whose ID_REF repeats a feature id (as gene symbols do
# in GSE136400: CXXC1 appears twice). Duplicate matrix rownames are tolerated
# but break ExpressionSet/AnnotatedDataFrame construction unless de-duplicated.
make_dup_series_matrix <- function() {
    txt <- c(
        "!Series_title\tdup-id fixture",
        "!Series_geo_accession\tGSE000098",
        "!Sample_title\tA\tB",
        "!Sample_geo_accession\tGSM000001\tGSM000002",
        "!Sample_platform_id\tGPL000001\tGPL000001",
        "!series_matrix_table_begin",
        "\"ID_REF\"\t\"GSM000001\"\t\"GSM000002\"",
        "\"CXXC1\"\t1.0\t2.0",
        "\"CXXC1\"\t3.0\t4.0",
        "\"CXXC4\"\t5.0\t6.0",
        "!series_matrix_table_end"
    )
    f <- tempfile(fileext = ".txt")
    writeLines(txt, f)
    f
}

test_that("parseGSEMatrix de-duplicates duplicate ID_REF feature IDs (#98)", {
    f <- make_dup_series_matrix()
    # getGPL = FALSE keeps this offline; the failure historically occurred while
    # building the ExpressionSet featureData, before any GPL fetch.
    eset <- GEOquery:::parseGSEMatrix(f, getGPL = FALSE)$eset

    expect_s4_class(eset, "ExpressionSet")
    expect_equal(Biobase::featureNames(eset), c("CXXC1", "CXXC1.1", "CXXC4"))
    expect_equal(nrow(Biobase::exprs(eset)), 3L)
    expect_equal(ncol(Biobase::exprs(eset)), 2L)
})

# --- #14: metadata-only ("quick"/"brief") SOFT files -----------------------

test_that("parseGSE returns a metadata-only GSE when there are no entities (#14)", {
    # The truncated SOFT from getGEOfile(amount = "quick") has only Series header
    # lines and no ^SAMPLE/^PLATFORM entities. This must parse, not error.
    txt <- c(
        "^SERIES = GSE000010",
        "!Series_title = quick-output fixture",
        "!Series_geo_accession = GSE000010",
        "!Series_summary = metadata only, no entities"
    )
    f <- tempfile(fileext = ".soft")
    writeLines(txt, f)

    gse <- suppressMessages(GEOquery:::parseGSE(f))
    expect_s4_class(gse, "GSE")
    expect_length(GSMList(gse), 0L)
    expect_length(GPLList(gse), 0L)
    expect_equal(Meta(gse)$title, "quick-output fixture")
})

# --- #80: getGSEDataTables single-table series -----------------------------

test_that(".parseGSEDataTableNodes returns a list of data.frames for one table (#80)", {
    # Mimic the acc.cgi form=xml Series document with a single <Data-Table>
    # (as GSE98638 has). sapply used to mis-simplify this into a matrix.
    doc <- xml2::read_xml(paste0(
        '<MINiML xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML">',
        '<Series>',
        '<Data-Table>',
        '<Column position="1"><Name>UniqueCell_ID</Name></Column>',
        '<Column position="2"><Name>Patient</Name></Column>',
        '<Internal-Data rows="2">',
        "NTC148-0322\tP0322\nNTC69-0407\tP0407",
        '</Internal-Data>',
        '</Data-Table>',
        '</Series>',
        '</MINiML>'
    ))
    nodes <- xml2::xml_find_all(doc, "//d1:Data-Table")
    res <- suppressWarnings(GEOquery:::.parseGSEDataTableNodes(nodes))

    expect_type(res, "list")
    expect_length(res, 1L)
    expect_s3_class(res[[1]], "data.frame")
    expect_equal(colnames(res[[1]]), c("UniqueCell_ID", "Patient"))
    expect_equal(nrow(res[[1]]), 2L)
})

test_that(".parseGSEDataTableNodes handles multiple tables (#80)", {
    doc <- xml2::read_xml(paste0(
        '<MINiML xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML">',
        '<Series>',
        '<Data-Table><Column position="1"><Name>a</Name></Column>',
        '<Internal-Data rows="2">x1\nx2</Internal-Data></Data-Table>',
        '<Data-Table><Column position="1"><Name>b</Name></Column>',
        '<Internal-Data rows="2">y1\ny2</Internal-Data></Data-Table>',
        '</Series></MINiML>'
    ))
    nodes <- xml2::xml_find_all(doc, "//d1:Data-Table")
    res <- suppressWarnings(GEOquery:::.parseGSEDataTableNodes(nodes))
    expect_type(res, "list")
    expect_length(res, 2L)
})

# --- #148: encoding override -----------------------------------------------

test_that(".geo_encoding() reflects the GEOquery.encoding option (#148)", {
    expect_equal(GEOquery:::.geo_encoding(), "unknown")
    old <- options(GEOquery.encoding = "Latin-1")
    on.exit(options(old), add = TRUE)
    expect_equal(GEOquery:::.geo_encoding(), "Latin-1")
})

test_that("parseGSEMatrix still parses under a non-default encoding option (#148)", {
    # A plain 2-feature series matrix (self-contained; no cross-file helpers).
    txt <- c(
        "!Series_title\tencoding fixture",
        "!Series_geo_accession\tGSE000148",
        "!Sample_title\tA\tB",
        "!Sample_geo_accession\tGSM000001\tGSM000002",
        "!Sample_platform_id\tGPL000001\tGPL000001",
        "!series_matrix_table_begin",
        "\"ID_REF\"\t\"GSM000001\"\t\"GSM000002\"",
        "\"probe_1\"\t1.5\t11.5",
        "\"probe_2\"\t2.5\t12.5",
        "!series_matrix_table_end"
    )
    f <- tempfile(fileext = ".txt")
    writeLines(txt, f)
    old <- options(GEOquery.encoding = "Latin-1")
    on.exit(options(old), add = TRUE)
    eset <- GEOquery:::parseGSEMatrix(f, getGPL = FALSE)$eset
    expect_s4_class(eset, "ExpressionSet")
    expect_equal(nrow(Biobase::exprs(eset)), 2L)
})

test_that("getGEO() rejects an invalid encoding argument (#148)", {
    expect_error(getGEO("GSE1", encoding = "bogus"), "should be one of")
})

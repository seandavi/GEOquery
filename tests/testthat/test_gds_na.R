# Offline regression test for #21: an NA in a GDS ID_REF column must not break
# GDS2eSet(). Build a minimal synthetic GDS object -- no network.

make_fake_gds <- function() {
    tbl <- data.frame(
        ID_REF = c("probe_1", "probe_2", NA),
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

test_that("GDS2eSet tolerates NA in ID_REF (#21)", {
    gds <- make_fake_gds()
    eset <- GDS2eSet(gds, getGPL = FALSE)

    expect_s4_class(eset, "ExpressionSet")
    expect_true(validObject(eset))
    expect_equal(nrow(Biobase::exprs(eset)), 3L)
    expect_equal(ncol(Biobase::exprs(eset)), 2L)
    # the NA ID_REF became a usable, non-missing feature name
    expect_false(any(is.na(Biobase::featureNames(eset))))
    expect_true("NA" %in% Biobase::featureNames(eset))
})

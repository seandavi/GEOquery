# Offline tests for the SummarizedExperiment return-type option (#168).

make_eset <- function(nrow = 3, ncol = 2) {
    m <- matrix(
        as.double(seq_len(nrow * ncol)), nrow = nrow,
        dimnames = list(paste0("g", seq_len(nrow)), paste0("s", seq_len(ncol)))
    )
    Biobase::ExpressionSet(assayData = m)
}

test_that("as_SummarizedExperiment coerces an ExpressionSet (#168)", {
    se <- as_SummarizedExperiment(make_eset(3, 2))
    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(dim(se), c(3L, 2L))
    expect_equal(unname(SummarizedExperiment::assay(se)[1, 1]), 1)
    expect_error(as_SummarizedExperiment("not an eset"), "ExpressionSet")
})

test_that(".applyReturnType coerces ExpressionSets and passes others through (#168)", {
    eset <- make_eset(2, 2)

    # SummarizedExperiment: list elements coerced
    out <- GEOquery:::.applyReturnType(list(a = eset), "SummarizedExperiment")
    expect_s4_class(out$a, "SummarizedExperiment")

    # ExpressionSet (default): unchanged
    out2 <- GEOquery:::.applyReturnType(list(a = eset), "ExpressionSet")
    expect_s4_class(out2$a, "ExpressionSet")

    # single ExpressionSet coerced
    expect_s4_class(
        GEOquery:::.applyReturnType(eset, "SummarizedExperiment"),
        "SummarizedExperiment"
    )

    # non-ExpressionSet objects (e.g. SOFT S4 results) pass through unchanged
    expect_equal(
        GEOquery:::.applyReturnType("soft-object", "SummarizedExperiment"),
        "soft-object"
    )
})

test_that("getGEO defaults to SummarizedExperiment (#168)", {
    # the default of the returnType argument is now SummarizedExperiment
    expect_equal(eval(formals(getGEO)$returnType)[1], "SummarizedExperiment")
})

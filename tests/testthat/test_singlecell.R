# Offline tests for single-cell file classification (#158, SC1). The pure
# classification logic is tested here; geoSingleCellManifest() itself hits the
# network and is covered by the integration tests.

test_that("single supplementary files are classified by format (#158)", {
    f <- GEOquery:::.classify_sc_file
    expect_equal(f("GSM123_matrix.mtx.gz")[["format"]], "10x_mtx")
    expect_equal(f("GSM123_matrix.mtx.gz")[["role"]], "matrix")
    expect_equal(f("GSM123_barcodes.tsv.gz")[["role"]], "barcodes")
    expect_equal(f("GSM123_features.tsv.gz")[["role"]], "features")
    expect_equal(f("GSM123_genes.tsv.gz")[["role"]], "features")
    expect_equal(f("study.h5ad")[["format"]], "h5ad")
    expect_equal(f("filtered_feature_bc_matrix.h5")[["format"]], "10x_h5")
    expect_equal(f("cells.loom")[["format"]], "loom")
    expect_equal(f("seurat.rds")[["format"]], "rds")
    expect_equal(f("GSE1_RAW.tar")[["format"]], "tar")
    expect_equal(f("readme.txt")[["format"]], "other")
})

test_that("the manifest groups by sample and reports format/role (#158)", {
    m <- GEOquery:::.classify_sc_files(
        c("GSM1_matrix.mtx.gz", "GSM1_barcodes.tsv.gz", "GSM1_features.tsv.gz", "GSM2_data.h5ad"),
        urls = paste0("https://x/", 1:4)
    )
    expect_equal(nrow(m), 4L)
    expect_equal(m$sample, c("GSM1", "GSM1", "GSM1", "GSM2"))
    expect_true(all(m$format[1:3] == "10x_mtx"))
    expect_equal(m$format[4], "h5ad")
    expect_equal(m$url[1], "https://x/1")
})

test_that(".sc_sample_id extracts GSM ids or NA", {
    expect_equal(GEOquery:::.sc_sample_id("GSM98765_matrix.mtx.gz"), "GSM98765")
    expect_true(is.na(GEOquery:::.sc_sample_id("combined_matrix.mtx.gz")))
})

test_that("geoSingleCellUnits groups by sample/format and assesses completeness (#158)", {
    manifest <- GEOquery:::.classify_sc_files(c(
        "GSM1_matrix.mtx.gz", "GSM1_barcodes.tsv.gz", "GSM1_features.tsv.gz", # complete 10x
        "GSM2_matrix.mtx.gz", "GSM2_barcodes.tsv.gz",                         # incomplete (no features)
        "GSM3_data.h5ad",                                                     # single-file complete
        "GSM4_notes.txt"                                                      # unsupported
    ))
    u <- geoSingleCellUnits(manifest)

    gsm1 <- u[u$sample == "GSM1", ]
    expect_equal(gsm1$status, "complete")
    expect_true(gsm1$loadable)
    expect_equal(gsm1$n_files, 3L)

    gsm2 <- u[u$sample == "GSM2", ]
    expect_match(gsm2$status, "incomplete")
    expect_match(gsm2$status, "features")
    expect_false(gsm2$loadable)

    expect_true(u[u$sample == "GSM3", "loadable"])
    expect_equal(u[u$sample == "GSM4", "status"], "unsupported")
    expect_false(u[u$sample == "GSM4", "loadable"])
})

test_that("geoSingleCellUnits handles an empty manifest", {
    empty <- GEOquery:::.classify_sc_files(character(0), character(0))
    u <- geoSingleCellUnits(empty)
    expect_equal(nrow(u), 0L)
    expect_true(all(c("sample", "format", "status", "loadable") %in% colnames(u)))
})

test_that("readGEOSingleCell rejects unsupported formats (#158)", {
    expect_error(readGEOSingleCell("x.loom", format = "loom"), "not supported")
    expect_error(readGEOSingleCell("x.rds", format = "rds"), "not supported")
})

test_that(".select_sc_units picks loadable units, one format per sample (#158)", {
    manifest <- GEOquery:::.classify_sc_files(c(
        "GSM1_data.h5ad",
        "GSM1_matrix.mtx.gz", "GSM1_barcodes.tsv.gz", "GSM1_features.tsv.gz",
        "GSM2_matrix.mtx.gz", "GSM2_barcodes.tsv.gz"   # incomplete
    ))
    units <- geoSingleCellUnits(manifest)
    sel <- GEOquery:::.select_sc_units(units)

    # GSM1 is loadable in two formats; h5ad wins by priority
    expect_equal(nrow(sel$load), 1L)
    expect_equal(sel$load$sample, "GSM1")
    expect_equal(sel$load$format, "h5ad")
    # GSM2 (incomplete) is among the skipped
    expect_true("GSM2" %in% sel$skip$sample)
})

test_that(".select_sc_units honors samples and format filters (#158)", {
    manifest <- GEOquery:::.classify_sc_files(c("GSM1_a.h5ad", "GSM2_b.h5ad"))
    units <- geoSingleCellUnits(manifest)

    expect_equal(GEOquery:::.select_sc_units(units, samples = "GSM2")$load$sample, "GSM2")
    expect_equal(nrow(GEOquery:::.select_sc_units(units, format = "10x_mtx")$load), 0L)
})

test_that("readGEOSingleCell reads a .rds SingleCellExperiment by class (#197)", {
    skip_if_not_installed("SingleCellExperiment")
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(as.double(1:6), nrow = 3))
    )
    f <- tempfile(fileext = ".rds")
    saveRDS(sce, f)

    out <- readGEOSingleCell(f)
    expect_s4_class(out, "SingleCellExperiment")
    expect_equal(dim(out), c(3L, 2L))
})

test_that("readGEOSingleCell rejects .rds with unsupported contents (#197)", {
    f <- tempfile(fileext = ".rds")
    saveRDS(matrix(1:4, 2), f)
    expect_error(readGEOSingleCell(f), "Unsupported .rds")
})

test_that("readGEOSingleCell still rejects loom; rds is now supported (#197)", {
    expect_error(readGEOSingleCell("x.loom", format = "loom"), "not supported")
    # the message no longer claims rds is unsupported
    err <- tryCatch(readGEOSingleCell("x.loom", format = "loom"), error = conditionMessage)
    expect_false(grepl("\\brds\\b.*not handled|Seurat .rds are not", err))
})

test_that(".as_sc_output passes SingleCellExperiment through (#196)", {
    skip_if_not_installed("SingleCellExperiment")
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(as.double(1:4), nrow = 2))
    )
    expect_identical(GEOquery:::.as_sc_output(sce, "SingleCellExperiment"), sce)
})

test_that(".as_sc_output coerces to Seurat when requested (#196)", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("Seurat")
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(as.double(1:6), nrow = 3,
            dimnames = list(c("g1", "g2", "g3"), c("c1", "c2"))))
    )
    obj <- GEOquery:::.as_sc_output(sce, "Seurat")
    expect_s4_class(obj, "Seurat")
})

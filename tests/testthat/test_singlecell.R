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

test_that("geoSingleCellManifest rejects non-GSE/GSM accessions", {
    expect_error(geoSingleCellManifest("GPL570"), "GSE or GSM")
    expect_error(geoSingleCellManifest("GDS507"), "GSE or GSM")
})

test_that(".sc_manifest_from_gsms returns an empty manifest for no input", {
    m <- GEOquery:::.sc_manifest_from_gsms(character(0))
    expect_equal(nrow(m), 0L)
    expect_true(all(c("fname", "sample", "format", "role", "url") %in% colnames(m)))
})

# A minimal SingleCellExperiment with the given feature (row) names. `rowdata`
# names the rowData columns to attach (each filled from the gene ids), mimicking
# 10x v2 (genes.tsv -> ID,Symbol) vs v3 (features.tsv -> ID,Symbol,Type).
.fake_sce <- function(genes, cells, rowdata = NULL) {
    m <- matrix(
        seq_len(length(genes) * length(cells)),
        nrow = length(genes),
        dimnames = list(genes, cells)
    )
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
    if (!is.null(rowdata)) {
        rd <- S4Vectors::DataFrame(
            lapply(stats::setNames(rowdata, rowdata), function(col) genes),
            row.names = genes
        )
        SummarizedExperiment::rowData(sce) <- rd
    }
    sce
}

test_that(".combine_sce binds samples that share all features (#190)", {
    skip_if_not_installed("SingleCellExperiment")
    a <- .fake_sce(c("g1", "g2", "g3"), c("c1", "c2"))
    b <- .fake_sce(c("g1", "g2", "g3"), c("c3", "c4", "c5"))
    out <- GEOquery:::.combine_sce(list(A = a, B = b))
    expect_equal(nrow(out), 3L)
    expect_equal(ncol(out), 5L)
    expect_equal(rownames(out), c("g1", "g2", "g3"))
})

test_that(".combine_sce restricts to common features when they differ (#190)", {
    skip_if_not_installed("SingleCellExperiment")
    a <- .fake_sce(c("g1", "g2", "g3"), c("c1", "c2"))
    b <- .fake_sce(c("g2", "g3", "g4"), c("c3", "c4"))
    expect_message(
        out <- GEOquery:::.combine_sce(list(A = a, B = b)),
        "Combining on 2 feature"
    )
    expect_equal(nrow(out), 2L)
    expect_equal(ncol(out), 4L)
    expect_setequal(rownames(out), c("g2", "g3"))
})

test_that(".combine_sce errors when samples share no features (#190)", {
    skip_if_not_installed("SingleCellExperiment")
    a <- .fake_sce(c("g1", "g2"), c("c1"))
    b <- .fake_sce(c("g3", "g4"), c("c2"))
    expect_error(
        GEOquery:::.combine_sce(list(A = a, B = b)),
        "no common features"
    )
})

test_that(".combine_sce binds samples with differing rowData columns (#190)", {
    skip_if_not_installed("SingleCellExperiment")
    # Same features, but 10x v2 (ID,Symbol) vs v3 (ID,Symbol,Type) rowData -- a
    # naive cbind() fails with "subscript contains invalid names". The combine
    # reduces to the shared rowData columns and binds.
    genes <- c("g1", "g2", "g3")
    a <- .fake_sce(genes, c("c1", "c2"), rowdata = c("ID", "Symbol"))
    b <- .fake_sce(genes, c("c3", "c4"), rowdata = c("ID", "Symbol", "Type"))
    out <- GEOquery:::.combine_sce(list(A = a, B = b))
    expect_equal(nrow(out), 3L)
    expect_equal(ncol(out), 4L)
    expect_equal(
        colnames(SummarizedExperiment::rowData(out)),
        c("ID", "Symbol")
    )
})

# ---- Integration tests (live network; GSE132771) --------------------------
# GSE132771 ships only a GSE..._RAW.tar at the series level; its per-sample 10x
# triplets live in each GSM suppl directory. Exercises the GSM-level fallback
# and the GSM-accession entry point.

test_that("geoSingleCellManifest(GSM) inventories a single sample (#190)", {
    skip_if_no_integration()
    m <- geoSingleCellManifest("GSM3891612")
    expect_equal(nrow(m), 3L)
    expect_true(all(m$sample == "GSM3891612"))
    expect_true(all(m$format == "10x_mtx"))
    expect_setequal(m$role, c("matrix", "barcodes", "features"))
})

test_that("geoSingleCellManifest(GSE) falls back to GSM-level files (#190)", {
    skip_if_no_integration()
    # Series level is only a _RAW.tar -> no loadable units -> GSM fallback.
    m <- geoSingleCellManifest("GSE132771")
    expect_gt(nrow(m), 3L)
    expect_true(all(m$format == "10x_mtx"))
    u <- geoSingleCellUnits(m)
    expect_true(all(u$loadable))
    expect_gt(sum(u$loadable), 1L)
})

test_that("geoSingleCellManifest(GSE, samples=) restricts without enumerating all (#190)", {
    skip_if_no_integration()
    m <- geoSingleCellManifest("GSE132771", samples = c("GSM3891612", "GSM3891613"))
    expect_setequal(unique(m$sample), c("GSM3891612", "GSM3891613"))
})

test_that("getGEOSingleCell reads a GSM into a SingleCellExperiment (#190)", {
    skip_if_no_integration()
    skip_if_not_installed("TENxIO")
    skip_if_not_installed("SingleCellExperiment")
    res <- getGEOSingleCell("GSM3891612", destdir = tempfile("sc_"))
    expect_type(res, "list")
    expect_named(res, "GSM3891612")
    expect_s4_class(res[[1]], "SingleCellExperiment")
    expect_gt(nrow(res[[1]]), 0L)
    expect_gt(ncol(res[[1]]), 0L)
})

test_that("getGEOSingleCell(GSE, samples=) reads selected samples (#190)", {
    skip_if_no_integration()
    skip_if_not_installed("TENxIO")
    skip_if_not_installed("SingleCellExperiment")
    res <- getGEOSingleCell("GSE132771",
        samples = c("GSM3891614", "GSM3891615"), destdir = tempfile("sc_"))
    expect_named(res, c("GSM3891614", "GSM3891615"))
    expect_s4_class(res[[1]], "SingleCellExperiment")
})

test_that("getGEOSingleCell(combine=TRUE) binds same-platform samples (#190)", {
    skip_if_no_integration()
    skip_if_not_installed("TENxIO")
    skip_if_not_installed("SingleCellExperiment")
    # GSM3891614/615 are both mouse (GPL21103) -> identical features.
    out <- getGEOSingleCell("GSE132771",
        samples = c("GSM3891614", "GSM3891615"), destdir = tempfile("sc_"),
        combine = TRUE)
    expect_s4_class(out, "SingleCellExperiment")
    expect_gt(ncol(out), 0L)
})

test_that("getGEOSingleCell(combine=TRUE) errors on incompatible platforms (#190)", {
    skip_if_no_integration()
    skip_if_not_installed("TENxIO")
    skip_if_not_installed("SingleCellExperiment")
    # GSM3891612 is mouse (GPL21103), GSM3891620 is human (GPL24676): no shared
    # features, so a combined object is impossible -- expect a clear error.
    expect_error(
        getGEOSingleCell("GSE132771",
            samples = c("GSM3891612", "GSM3891620"), destdir = tempfile("sc_"),
            combine = TRUE),
        "no common features"
    )
})

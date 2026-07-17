# Offline tests for per-sample GEO metadata -> colData (#210). The metadata
# *fetch* helpers (.gse_sample_metadata / .gsm_sample_metadata) hit the network
# and are exercised by the integration suite; here we test the pure logic and
# the SCE-manipulating helpers with constructed objects.

test_that(".sample_meta_cols selects characteristics + title, not technical fields (#210)", {
    cn <- c(
        "title", "geo_accession", "source_name_ch1", "organism_ch1",
        "characteristics_ch1", "characteristics_ch1.1", "molecule_ch1",
        "data_processing", "contact_name", "platform_id",
        "age.ch1", "Sex.ch1", "genotype.ch1", "tissue.ch1", "treatment.ch2"
    )
    got <- GEOquery:::.sample_meta_cols(cn)
    expect_setequal(
        got,
        c("title", "source_name_ch1", "age.ch1", "Sex.ch1", "genotype.ch1",
            "tissue.ch1", "treatment.ch2")
    )
    # raw/technical fields are excluded
    expect_false(any(c("characteristics_ch1", "organism_ch1", "molecule_ch1",
        "data_processing", "contact_name", "platform_id") %in% got))
})

test_that(".split_characteristics parses key: value strings (#210)", {
    kv <- GEOquery:::.split_characteristics(c("age: 2 month", "Sex: Male", "", NA, "no colon"))
    expect_equal(kv[["age"]], "2 month")
    expect_equal(kv[["Sex"]], "Male")
    expect_length(kv, 2L)
    # duplicate keys are made unique
    kv2 <- GEOquery:::.split_characteristics(c("rep: 1", "rep: 2"))
    expect_setequal(names(kv2), c("rep", "rep.1"))
    expect_length(GEOquery:::.split_characteristics(character(0)), 0L)
})

test_that(".attach_sample_meta broadcasts a metadata row across cells, prefixed (#210)", {
    skip_if_not_installed("SingleCellExperiment")
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(as.double(1:6), nrow = 3,
            dimnames = list(c("g1", "g2", "g3"), c("c1", "c2"))))
    )
    meta <- data.frame(title = "SampleA", age.ch1 = "2 month",
        stringsAsFactors = FALSE, check.names = FALSE)
    out <- GEOquery:::.attach_sample_meta(sce, meta)
    cd <- SummarizedExperiment::colData(out)
    expect_true(all(c("sample.title", "sample.age.ch1") %in% colnames(cd)))
    expect_equal(as.character(cd$sample.title), c("SampleA", "SampleA"))
    expect_equal(as.character(cd[["sample.age.ch1"]]), c("2 month", "2 month"))
    # a NULL/empty metadata row leaves the object untouched
    expect_identical(GEOquery:::.attach_sample_meta(sce, NULL), sce)
})

test_that(".combine_sce reconciles differing colData columns (#210)", {
    skip_if_not_installed("SingleCellExperiment")
    mk <- function(cells, meta_col, val) {
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = matrix(as.double(seq_len(3 * length(cells))),
                nrow = 3, dimnames = list(c("g1", "g2", "g3"), cells)))
        )
        SummarizedExperiment::colData(sce)[[meta_col]] <- rep(val, length(cells))
        sce
    }
    a <- mk(c("a1", "a2"), "sample.age.ch1", "2 month")   # has age, not sex
    b <- mk(c("b1", "b2"), "sample.Sex.ch1", "Male")      # has sex, not age
    combined <- GEOquery:::.combine_sce(list(A = a, B = b))

    expect_equal(ncol(combined), 4L)
    cd <- SummarizedExperiment::colData(combined)
    expect_true(all(c("sample.age.ch1", "sample.Sex.ch1") %in% colnames(cd)))
    # A's cells have age set, sex NA; B's cells the reverse
    expect_equal(as.character(cd[c("a1", "a2"), "sample.age.ch1"]), c("2 month", "2 month"))
    expect_true(all(is.na(cd[c("a1", "a2"), "sample.Sex.ch1"])))
    expect_equal(as.character(cd[c("b1", "b2"), "sample.Sex.ch1"]), c("Male", "Male"))
    expect_true(all(is.na(cd[c("b1", "b2"), "sample.age.ch1"])))
})

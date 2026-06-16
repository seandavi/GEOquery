# GEOquery 2.81.23 (2026-06-16)

## New features

- Seurat interoperability for single-cell data (optional; `Seurat` in `Suggests`). `readGEOSingleCell()` now reads `.rds` supplementary files containing a Seurat or `SingleCellExperiment` object (detected by class; Seurat coerced to `SingleCellExperiment`), so `.rds` is now a loadable single-cell format. Both `readGEOSingleCell()` and `getGEOSingleCell()` gain `as = "Seurat"` to return Seurat objects (coerced at the output boundary). `SingleCellExperiment` remains the internal representation; see ADR-0006 (#195, #196, #197).

# GEOquery 2.81.22 (2026-06-14)

## New features

- `geoSingleCellManifest()` and `getGEOSingleCell()` now handle the common case where a Series ships only a `GSE..._RAW.tar` at the series level and the per-sample files live in each sample's own GSM suppl directory (e.g. GSE132771). When the series level has no loadable single-cell units, the manifest falls back to enumerating the Series' samples (via `getGEO()`) and inventorying each GSM suppl directory. Both functions also accept a GSM accession directly (`getGEOSingleCell("GSM3891612")`), and `geoSingleCellManifest()` gains a `samples` argument to restrict to specific GSMs without enumerating the whole Series. Unit files are downloaded by URL, so the readers work whether the data lives at the series or sample level (#190).
- `getGEOSingleCell()` gains a `by` argument — one of `"sample"` (default; a named list of one `SingleCellExperiment` per sample), `"platform"` (a named list keyed by platform/GPL, each platform's samples combined), or `"all"` (a single combined object) — replacing the previous `combine` flag. A GEO Series can span multiple platforms (e.g. GSE132771 mixes mouse and human), and the platform is the natural feature-compatibility boundary, so `by = "platform"` is the right way to combine a multi-platform study; the return shape is determined by the argument, not the data. `geoSingleCellManifest()` now reports a `platform` (GPL) column per sample, and `geoSingleCellUnits()` treats each whole-study single file as its own unit so a study's several `.h5ad` files are no longer merged into one bogus unit. Gzipped single-file formats (`.h5ad.gz`, `.h5.gz`) are recognized and transparently decompressed before reading. loom is reported but flagged not loadable, no built-in reader (`.rds`/Seurat support follows in 2.81.23) (#190).

## Bug Fixes

- Combining single-cell samples (`getGEOSingleCell(by = "platform"|"all")`) no longer fails with a cryptic `cbind` error (`'mcols' ... do not match` or `subscript contains invalid names`) when a Series' samples have heterogeneous feature annotation — common in single-cell studies (e.g. GSE132771 mixes 10x CellRanger v2 `genes.tsv` with v3 `features.tsv`, giving different rowData columns). Samples are restricted to their shared features and given one canonical rowData before binding; when they share no features (so a single combined object is impossible) a clear, actionable error is raised (#190).

# GEOquery 2.81.21 (2026-06-13)

## Breaking changes

- **`getGEO()` now returns `SummarizedExperiment` objects by default** for GSE Series Matrix records (previously `ExpressionSet`). Update downstream code from `exprs()`/`pData()`/`fData()` to `assay()`/`colData()`/`rowData()`, or pass `returnType = "ExpressionSet"` to keep the old behavior. Existing results can also be converted with `as_SummarizedExperiment()`. SOFT-format results (GDS/GPL/GSM/GSE S4 objects) are unaffected. See ADR-0002 and ADR-0005 (#168).

## New features

- Optional persistent download cache backed by **BiocFileCache**. Set `options(GEOquery.cache = TRUE)` to have downloads keyed on their URL and reused across sessions (location defaults to `tools::R_user_dir("GEOquery", "cache")`, overridable via `options(GEOquery.cache.path = ...)`). New `geoCache()` and `clearGEOCache()` expose and clear it. Off by default for now, preserving the historical `destdir` behavior (#171).
- New `readGEOSingleCell()` and `getGEOSingleCell()` read GEO single-cell supplementary data into `SingleCellExperiment` objects: 10x Matrix Market and 10x HDF5 via **TENxIO**, AnnData `.h5ad` via **anndataR** (optional `Suggests`). `getGEOSingleCell()` returns a named list of per-sample objects (combine with care) and reports which units it loads and skips. loom, files inside `_RAW.tar`, and idiosyncratic layouts are intentionally out of scope — use `geoSingleCellManifest()` + `readGEOSingleCell()` for those (#158, #190).
- New `geoSingleCellManifest()` inventories a GSE's supplementary files and classifies them by single-cell format (10x Matrix Market triplet, 10x HDF5, AnnData h5ad, loom, Seurat rds, tar), grouping by GSM sample — so you can see what a single-cell study contains before downloading. `geoSingleCellUnits()` collapses the manifest into loadable units (per sample + format) and flags completeness (e.g. an incomplete 10x triplet). Steps toward single-cell readers (ADR-0004) (#158, #188, #189).
- `getGEO()` gains a `returnType` argument. With `returnType = "SummarizedExperiment"`, GSE Series Matrix results are returned as `SummarizedExperiment` objects instead of `ExpressionSet`. The default remains `"ExpressionSet"` for now (with a one-time notice) and will switch to `"SummarizedExperiment"` in a future release. A new exported `as_SummarizedExperiment()` coerces an existing `ExpressionSet` result without re-downloading. See ADR-0002 (#168).
- Downloads now stream to disk instead of buffering the entire response in memory, retry on transient HTTP errors, and honor a configurable `GEOquery.download.timeout` option (default 300 seconds) — replacing the previous enforced 120-second floor that ignored lower user timeouts. Failures raise a typed `geoquery_download_error` carrying the URL and HTTP status. `getDirListing()` now uses the same httr2 layer (#147, #173).
- GEOquery now raises typed error conditions — `geoquery_error` and subclasses (`geoquery_private_accession`, `geoquery_download_error`, `geoquery_parse_error`, `geoquery_bad_accession`) — so failures can be handled programmatically with `tryCatch()` (#170, #184, #186).
- `getGEOSuppFiles()` gains a `quiet` argument (defaulting to the `GEOquery.quiet` option, or `FALSE`) to suppress informational messages such as "No supplemental files found" and "Using locally cached version" (#68, #182).

## Documentation

- The S4 class and accessor documentation is filled in: the `GEOData` accessors (`Meta`, `Table`, `Columns`, `dataTable`, `Accession`, `GSMList`, `GPLList`) now have real descriptions, return values, and examples, and the class pages no longer imply constructing objects with `new()` — they are returned by `getGEO()` (#103, #192).
- Documentation is reorganized into narrative pkgdown **articles** — *Understanding GEO data formats*, *RNA-seq quantifications*, *Single-cell data from GEO*, and *From GEO to downstream analysis* — that cover the *why* (entity types, file formats) and downstream workflows with links to other Bioconductor packages. The package vignette is now a concise quick-start that indexes them; the articles render on the pkgdown site and are excluded from `R CMD check` (#156, #191).

- The package `DESCRIPTION` and `biocViews` now describe GEOquery's actual scope (microarray, RNA-seq, and single-cell; GEO Series Matrix files parsed to `ExpressionSet` by default) instead of microarray-only (#71, #181).

## Bug Fixes

- Supplemental-file URLs are now built with a small `url_join()` helper instead of `file.path()`, which mangled `https://` into `https:/` and produced double slashes. Affects `getGEOSuppFiles(fetch_files = FALSE)` and `getGEOSeriesFileListing()` (#131, #178).
- `GDS2eSet()` no longer fails when a GDS has an `NA` (or empty) value in its `ID_REF` column (e.g. GDS3666). Such values are replaced with a usable feature name instead of producing "row names contain missing values" (#21, #177).
- `getGEO()` now fails with a clear message when an accession is private, embargoed, or not yet public (NCBI returns an HTML page) instead of mis-parsing it or, in older versions, looping. `findFirstEntity()` is also hardened against a multi-line edge case that could error and against unbounded reads (#58, #176).
- `getGEO(parseCharacteristics = FALSE)` now actually skips characteristics parsing. The flag was accepted at the top level but dropped before reaching `parseGSEMatrix()`; it is now threaded through `getAndParseGSEMatrices()` and `parseGEO()` (#60, #175).
- Fixed error when parsing GSE matrix files with malformed or empty lines between sample metadata (e.g., GSE425). Sample lines are now extracted directly using pattern matching to avoid issues with irregular file formatting (#162).

# GEOquery 2.75.0 (2024-10-01)

## New Features

- RNAseq data support for GEOquery. Now you can use RNASeq quantification data prepared by NCBI.
- Basic search in GEO database. Now you can search for datasets in GEO database using GEOquery.
- browseGEO() function to open a web browser with a GEO accession.

## Bug Fixes or Improvements

Not an exhaustive list, but some highlights:

- Using httr2 instead of curl for better control over HTTP requests.
- Removed dead gunzip code.


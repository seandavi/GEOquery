# Changelog

## GEOquery (development version)

### New features

- New
  [`readGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md)
  and
  [`getGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/getGEOSingleCell.md)
  read GEO single-cell supplementary data into `SingleCellExperiment`
  objects: 10x Matrix Market and 10x HDF5 via **TENxIO**, AnnData
  `.h5ad` via **anndataR** (optional `Suggests`).
  [`getGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/getGEOSingleCell.md)
  returns a named list of per-sample objects (combine with care) and
  reports which units it loads and skips. loom, Seurat `.rds`, files
  inside `_RAW.tar`, and idiosyncratic layouts are intentionally out of
  scope — use
  [`geoSingleCellManifest()`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md) +
  [`readGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md)
  for those ([\#158](https://github.com/seandavi/GEOquery/issues/158),
  [\#190](https://github.com/seandavi/GEOquery/issues/190)).
- New
  [`geoSingleCellManifest()`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md)
  inventories a GSE’s supplementary files and classifies them by
  single-cell format (10x Matrix Market triplet, 10x HDF5, AnnData h5ad,
  loom, Seurat rds, tar), grouping by GSM sample — so you can see what a
  single-cell study contains before downloading.
  [`geoSingleCellUnits()`](http://seandavi.github.io/GEOquery/reference/geoSingleCellUnits.md)
  collapses the manifest into loadable units (per sample + format) and
  flags completeness (e.g. an incomplete 10x triplet). Steps toward
  single-cell readers (ADR-0004)
  ([\#158](https://github.com/seandavi/GEOquery/issues/158),
  [\#188](https://github.com/seandavi/GEOquery/issues/188),
  [\#189](https://github.com/seandavi/GEOquery/issues/189)).
- [`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
  gains a `returnType` argument. With
  `returnType = "SummarizedExperiment"`, GSE Series Matrix results are
  returned as `SummarizedExperiment` objects instead of `ExpressionSet`.
  The default remains `"ExpressionSet"` for now (with a one-time notice)
  and will switch to `"SummarizedExperiment"` in a future release. A new
  exported
  [`as_SummarizedExperiment()`](http://seandavi.github.io/GEOquery/reference/as_SummarizedExperiment.md)
  coerces an existing `ExpressionSet` result without re-downloading. See
  ADR-0002 ([\#168](https://github.com/seandavi/GEOquery/issues/168)).
- Downloads now stream to disk instead of buffering the entire response
  in memory, retry on transient HTTP errors, and honor a configurable
  `GEOquery.download.timeout` option (default 300 seconds) — replacing
  the previous enforced 120-second floor that ignored lower user
  timeouts. Failures raise a typed `geoquery_download_error` carrying
  the URL and HTTP status.
  [`getDirListing()`](http://seandavi.github.io/GEOquery/reference/getDirListing.md)
  now uses the same httr2 layer
  ([\#147](https://github.com/seandavi/GEOquery/issues/147),
  [\#173](https://github.com/seandavi/GEOquery/issues/173)).
- GEOquery now raises typed error conditions — `geoquery_error` and
  subclasses (`geoquery_private_accession`, `geoquery_download_error`,
  `geoquery_parse_error`, `geoquery_bad_accession`) — so failures can be
  handled programmatically with
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html)
  ([\#170](https://github.com/seandavi/GEOquery/issues/170),
  [\#184](https://github.com/seandavi/GEOquery/issues/184),
  [\#186](https://github.com/seandavi/GEOquery/issues/186)).
- [`getGEOSuppFiles()`](http://seandavi.github.io/GEOquery/reference/getGEOSuppFiles.md)
  gains a `quiet` argument (defaulting to the `GEOquery.quiet` option,
  or `FALSE`) to suppress informational messages such as “No
  supplemental files found” and “Using locally cached version”
  ([\#68](https://github.com/seandavi/GEOquery/issues/68),
  [\#182](https://github.com/seandavi/GEOquery/issues/182)).

### Documentation

- The S4 class and accessor documentation is filled in: the `GEOData`
  accessors (`Meta`, `Table`, `Columns`, `dataTable`, `Accession`,
  `GSMList`, `GPLList`) now have real descriptions, return values, and
  examples, and the class pages no longer imply constructing objects
  with [`new()`](https://rdrr.io/r/methods/new.html) — they are returned
  by
  [`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
  ([\#103](https://github.com/seandavi/GEOquery/issues/103),
  [\#192](https://github.com/seandavi/GEOquery/issues/192)).

- Documentation is reorganized into narrative pkgdown **articles** —
  *Understanding GEO data formats*, *RNA-seq quantifications*,
  *Single-cell data from GEO*, and *From GEO to downstream analysis* —
  that cover the *why* (entity types, file formats) and downstream
  workflows with links to other Bioconductor packages. The package
  vignette is now a concise quick-start that indexes them; the articles
  render on the pkgdown site and are excluded from `R CMD check`
  ([\#156](https://github.com/seandavi/GEOquery/issues/156),
  [\#191](https://github.com/seandavi/GEOquery/issues/191)).

- The package `DESCRIPTION` and `biocViews` now describe GEOquery’s
  actual scope (microarray, RNA-seq, and single-cell; GEO Series Matrix
  files parsed to `ExpressionSet` by default) instead of microarray-only
  ([\#71](https://github.com/seandavi/GEOquery/issues/71),
  [\#181](https://github.com/seandavi/GEOquery/issues/181)).

### Bug Fixes

- Supplemental-file URLs are now built with a small `url_join()` helper
  instead of [`file.path()`](https://rdrr.io/r/base/file.path.html),
  which mangled `https://` into `https:/` and produced double slashes.
  Affects `getGEOSuppFiles(fetch_files = FALSE)` and
  [`getGEOSeriesFileListing()`](http://seandavi.github.io/GEOquery/reference/getGEOSeriesFileListing.md)
  ([\#131](https://github.com/seandavi/GEOquery/issues/131),
  [\#178](https://github.com/seandavi/GEOquery/issues/178)).
- [`GDS2eSet()`](http://seandavi.github.io/GEOquery/reference/coercion.md)
  no longer fails when a GDS has an `NA` (or empty) value in its
  `ID_REF` column (e.g. GDS3666). Such values are replaced with a usable
  feature name instead of producing “row names contain missing values”
  ([\#21](https://github.com/seandavi/GEOquery/issues/21),
  [\#177](https://github.com/seandavi/GEOquery/issues/177)).
- [`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
  now fails with a clear message when an accession is private,
  embargoed, or not yet public (NCBI returns an HTML page) instead of
  mis-parsing it or, in older versions, looping. `findFirstEntity()` is
  also hardened against a multi-line edge case that could error and
  against unbounded reads
  ([\#58](https://github.com/seandavi/GEOquery/issues/58),
  [\#176](https://github.com/seandavi/GEOquery/issues/176)).
- `getGEO(parseCharacteristics = FALSE)` now actually skips
  characteristics parsing. The flag was accepted at the top level but
  dropped before reaching
  [`parseGSEMatrix()`](http://seandavi.github.io/GEOquery/reference/parseGSEMatrix.md);
  it is now threaded through `getAndParseGSEMatrices()` and
  [`parseGEO()`](http://seandavi.github.io/GEOquery/reference/parseGEO.md)
  ([\#60](https://github.com/seandavi/GEOquery/issues/60),
  [\#175](https://github.com/seandavi/GEOquery/issues/175)).
- Fixed error when parsing GSE matrix files with malformed or empty
  lines between sample metadata (e.g., GSE425). Sample lines are now
  extracted directly using pattern matching to avoid issues with
  irregular file formatting.

## GEOquery 2.75.0 (2024-10-01)

### New Features

- RNAseq data support for GEOquery. Now you can use RNASeq
  quantification data prepared by NCBI.
- Basic search in GEO database. Now you can search for datasets in GEO
  database using GEOquery.
- browseGEO() function to open a web browser with a GEO accession.

### Bug Fixes or Improvements

Not an exhaustive list, but some highlights:

- Using httr2 instead of curl for better control over HTTP requests.
- Removed dead gunzip code.

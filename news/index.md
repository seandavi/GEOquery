# Changelog

## GEOquery (development version)

### New features

- [`getGEOSuppFiles()`](http://seandavi.github.io/GEOquery/reference/getGEOSuppFiles.md)
  gains a `quiet` argument (defaulting to the `GEOquery.quiet` option,
  or `FALSE`) to suppress informational messages such as “No
  supplemental files found” and “Using locally cached version”
  ([\#68](https://github.com/seandavi/GEOquery/issues/68),
  [\#182](https://github.com/seandavi/GEOquery/issues/182)).

### Documentation

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

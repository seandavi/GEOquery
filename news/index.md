# Changelog

## GEOquery 2.99.1 (unreleased)

### Bug Fixes

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

## GEOquery 2.99.0 (2024-10-01)

### New Features

- RNAseq data support for GEOquery. Now you can use RNASeq
  quantification data prepared by NCBI.
- Basic search in GEO database. Now you can search for datasets in GEO
  database using GEOquery.
- browseGEO() function to open a web browser with a GEO accession.

### Breaking changes

- [`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
  now returns a list of SummarizedExperiment objects. This is a breaking
  change from previous versions of GEOquery. If you are using GEOquery
  in a script, you will need to update your code to reflect this change.

### Bug Fixes or Improvements

Not an exhaustive list, but some highlights:

- Using httr2 instead of curl for better control over HTTP requests.
- Removed dead gunzip code.

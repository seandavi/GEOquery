# GEOquery 2.99.1 (unreleased)

## Bug Fixes

- `getGEO()` now fails with a clear message when an accession is private, embargoed, or not yet public (NCBI returns an HTML page) instead of mis-parsing it or, in older versions, looping. `findFirstEntity()` is also hardened against a multi-line edge case that could error and against unbounded reads (#58, #176).
- `getGEO(parseCharacteristics = FALSE)` now actually skips characteristics parsing. The flag was accepted at the top level but dropped before reaching `parseGSEMatrix()`; it is now threaded through `getAndParseGSEMatrices()` and `parseGEO()` (#60, #175).
- Fixed error when parsing GSE matrix files with malformed or empty lines between sample metadata (e.g., GSE425). Sample lines are now extracted directly using pattern matching to avoid issues with irregular file formatting.

# GEOquery 2.99.0 (2024-10-01)

## New Features

- RNAseq data support for GEOquery. Now you can use RNASeq quantification data prepared by NCBI.
- Basic search in GEO database. Now you can search for datasets in GEO database using GEOquery.
- browseGEO() function to open a web browser with a GEO accession.

## Breaking changes

- `getGEO()` now returns a list of SummarizedExperiment objects. This is a breaking change from previous versions of GEOquery. If you are using GEOquery in a script, you will need to update your code to reflect this change.

## Bug Fixes or Improvements

Not an exhaustive list, but some highlights:

- Using httr2 instead of curl for better control over HTTP requests.
- Removed dead gunzip code.


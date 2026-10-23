# GEOquery 2.80.1

## Bug Fixes

- `getRNASeqData()` failed with "invalid 'row.names' length" because NCBI now
  serves the shared gene annotation table behind a reCAPTCHA page. The
  annotation is now read from the GEO FTP mirror
  (`geo/rnaseq_counts/annot/`) when available, and a non-table response
  gives a clear error instead (#235).

# GEOquery 2.99.1 (unreleased)

## Bug Fixes

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


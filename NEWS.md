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


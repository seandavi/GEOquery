# GEOquery 2.80.1

## Bug fixes

- The `GEOquery.qmd` vignette is now precompiled, so `R CMD build` makes no
  NCBI calls. The live source is `vignettes/GEOquery.qmd.orig`; regenerate the
  shipped `vignettes/GEOquery.qmd` with `Rscript dev/precompute-vignettes.R`.
  This fixes the Bioconductor 3.23 build errors in which the vignette died on
  transient NCBI failures: an E-utilities HTTP 500 inside `searchGEO()` on
  `nebbiolo1`, and an FTP `SSL_ERROR_SYSCALL` while `getGEO("GSE2553")` was
  downloading on `kunpeng2`. This is a backport of the precompiled-vignette
  approach from devel (2.81.28).
- The `searchFieldsGEO()`, `hasRNASeqQuantifications()`,
  `getRNASeqQuantGenomeInfo()`, and `getRNASeqData()` examples are now in
  `\donttest{}`, so `R CMD check` does not depend on NCBI E-utilities or the
  GEO download service being healthy. Devel already had the
  `searchFieldsGEO()` failure on a builder (EInfo returned HTML instead of
  XML).

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


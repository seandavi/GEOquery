# Getting started with GEOquery

GEOquery is the bridge between the NCBI Gene Expression Omnibus (GEO)
and Bioconductor: it downloads and parses GEO records into Bioconductor
objects. This page is a short quick-start and an index; the in-depth,
narrative documentation lives in the **articles** listed below.

## Install

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("GEOquery")
```

## Quick start

``` r

library(GEOquery)

# A GSE via the fast Series Matrix path -> a list of ExpressionSet,
# one per platform.
gse <- getGEO("GSE2553")
eset <- gse[[1]]
exprs(eset)    # expression matrix
pData(eset)    # sample metadata
fData(eset)    # feature annotation

# Other entity types parse to GEOquery's S4 classes:
getGEO("GSM11805")   # a sample
getGEO("GPL96")      # a platform
getGEO("GDS507")     # a curated dataset

# See what supplementary files a study has, without downloading:
getGEOSuppFiles("GSE63137", fetch_files = FALSE)
```

## In-depth articles

The articles go beyond the *how* to the *why* — the structure of GEO,
the file formats, and how a GEOquery object connects to downstream
Bioconductor workflows:

- [**Understanding GEO data
  formats**](https://seandavi.github.io/GEOquery/articles/geo-data-formats.html)
  — the four entity types, SOFT vs. Series Matrix, why
  [`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
  returns different classes, and `ExpressionSet`
  vs. `SummarizedExperiment`.
- [**RNA-seq quantifications from
  GEO**](https://seandavi.github.io/GEOquery/articles/rnaseq.html) —
  NCBI’s uniformly-computed counts and how to retrieve them.
- [**Single-cell data from
  GEO**](https://seandavi.github.io/GEOquery/articles/single-cell.html)
  — why single-cell data lives in supplementary files, and the inspect →
  decide → load workflow into a `SingleCellExperiment`.
- [**From GEO to downstream
  analysis**](https://seandavi.github.io/GEOquery/articles/downstream-analysis.html)
  — taking a GEOquery object into limma / DESeq2 / edgeR / the
  single-cell ecosystem, with links to the relevant packages.

## Getting help

- Usage questions: the [Bioconductor support
  site](https://support.bioconductor.org/), tagged `geoquery`.
- Bugs and feature requests: the [issue
  tracker](https://github.com/seandavi/GEOquery/issues) — please include
  a GEO accession and
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).

# GEOquery

> The bridge between the NCBI Gene Expression Omnibus (GEO) and Bioconductor.

<!-- badges: start -->
[![R-CMD-check](https://github.com/seandavi/GEOquery/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seandavi/GEOquery/actions/workflows/R-CMD-check.yaml)
[![Bioc release](https://bioconductor.org/shields/build/release/bioc/GEOquery.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/GEOquery)
[![Bioc devel](https://bioconductor.org/shields/build/devel/bioc/GEOquery.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/GEOquery)
[![Downloads](https://bioconductor.org/shields/downloads/release/GEOquery.svg)](https://bioconductor.org/packages/GEOquery)
[![Years in Bioc](https://bioconductor.org/shields/years-in-bioc/GEOquery.svg)](https://bioconductor.org/packages/GEOquery)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

GEOquery downloads and parses data from the NCBI [Gene Expression
Omnibus](https://www.ncbi.nlm.nih.gov/geo/) — a public repository of
high-throughput functional genomics data — into Bioconductor objects, so you
can go from a GEO accession to an analysis-ready object in one call.

## Capabilities

- **Series, Samples, Platforms, DataSets.** Parse any GEO entity (`GSE`, `GSM`,
  `GPL`, `GDS`) from either the compact Series Matrix or the full SOFT format.
- **Modern object model.** GSE Series Matrix records return
  `SummarizedExperiment` objects by default (or `ExpressionSet` via
  `returnType = "ExpressionSet"`).
- **RNA-seq.** Retrieve NCBI's uniformly-computed RNA-seq quantifications with
  `getRNASeqData()`.
- **Single-cell.** Inventory, group, and load single-cell supplementary data
  (10x Matrix Market, 10x HDF5, AnnData `.h5ad`, Seurat `.rds`) into
  `SingleCellExperiment` (or `Seurat`) objects.
- **Supplementary files.** List and download any attached files with
  `getGEOSuppFiles()`.
- **Search.** Query GEO programmatically with `searchGEO()`.
- **Robust downloads.** Streaming downloads with retries, an optional persistent
  `BiocFileCache` cache, and typed error conditions for `tryCatch()`.

## Installation

```r
# from Bioconductor (recommended)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("GEOquery")

# or the development version from GitHub
BiocManager::install("seandavi/GEOquery")
```

## Quick start

```r
library(GEOquery)

# A GSE via the fast Series Matrix path -> a list of SummarizedExperiment,
# one per platform.
gse <- getGEO("GSE2553")
se <- gse[[1]]
assay(se)      # expression matrix
colData(se)    # sample metadata
rowData(se)    # feature annotation

# Other entity types parse to GEOquery's S4 classes:
getGEO("GSM11805")   # a sample
getGEO("GPL96")      # a platform
getGEO("GDS507")     # a curated dataset

# See what supplementary files a study has, without downloading:
getGEOSuppFiles("GSE63137", fetch_files = FALSE)
```

## Documentation

The package vignette is a quick-start; the in-depth, narrative articles cover
the *why* and the downstream workflows:

- [Getting started](http://seandavi.github.io/GEOquery/)
- [Understanding GEO data formats](http://seandavi.github.io/GEOquery/articles/geo-data-formats.html)
- [RNA-seq quantifications](http://seandavi.github.io/GEOquery/articles/rnaseq.html)
- [Single-cell data from GEO](http://seandavi.github.io/GEOquery/articles/single-cell.html)
- [From GEO to downstream analysis](http://seandavi.github.io/GEOquery/articles/downstream-analysis.html)

Bioconductor landing pages:
[release](https://bioconductor.org/packages/release/bioc/html/GEOquery.html) ·
[devel](https://bioconductor.org/packages/devel/bioc/html/GEOquery.html)

## Getting help

- **Usage questions:** the [Bioconductor support
  site](https://support.bioconductor.org/), tagged `geoquery`.
- **Bugs and feature requests:** the [issue
  tracker](https://github.com/seandavi/GEOquery/issues) — please include a GEO
  accession and `sessionInfo()`.

## Contributing

Contributions are welcome as [pull
requests](https://github.com/seandavi/GEOquery/pulls) or [issues](https://github.com/seandavi/GEOquery/issues).
See [CONTRIBUTING.md](CONTRIBUTING.md) for the development workflow, and follow
the [Bioconductor coding standards](https://contributions.bioconductor.org/r-code.html)
where possible.

## Citation

If you use GEOquery, please cite:

> Davis S, Meltzer PS. *GEOquery: a bridge between the Gene Expression Omnibus
> (GEO) and BioConductor.* Bioinformatics. 2007;23(14):1846–1847.

```r
citation("GEOquery")
```

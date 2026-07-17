# Inventory the single-cell supplementary files of a GEO Series or Sample

Lists the supplementary files attached to a GSE (or a single GSM) and
classifies each by single-cell format (10x Matrix Market triplet, 10x
HDF5, AnnData h5ad, loom, Seurat rds, tar archive, or other), extracting
the GSM sample id where present. This lets you see what a single-cell
study contains – and how 10x triplets group by sample – before
downloading potentially many gigabytes.

## Usage

``` r
geoSingleCellManifest(GEO, samples = NULL)
```

## Arguments

- GEO:

  A GEO Series (`"GSE..."`) or Sample (`"GSM..."`) accession, e.g.
  "GSE132771" or "GSM3891612".

- samples:

  Optional character vector of GSM ids. For a GSE, restricts the
  inventory to these samples; when the series level has no loadable
  units this also avoids enumerating the whole series. Ignored when
  `GEO` is a GSM.

## Value

A data.frame with columns `fname`, `sample` (GSM id or NA), `format`,
`role`, and `url`. Zero rows if nothing is found.

## Details

For a GSE, the series-level suppl directory is inventoried first. Many
single-cell studies ship only a `GSE..._RAW.tar` there, with the
loadable per-sample files (10x triplets, h5, h5ad) living in each
sample's own GSM suppl directory. When the series level yields no
loadable units, the manifest falls back to enumerating the series'
samples (via
[`getGEO`](http://seandavi.github.io/GEOquery/reference/getGEO.md)) and
inventorying each GSM suppl directory. Pass a GSM accession to inventory
just that one sample.

No files are downloaded. The result feeds the single-cell readers (see
ADR-0004); reading itself uses Bioconductor importers (TENxIO, anndataR)
that are optional dependencies.

## See also

[`getGEOSuppFiles`](http://seandavi.github.io/GEOquery/reference/getGEOSuppFiles.md),
[`getGEOSingleCell`](http://seandavi.github.io/GEOquery/reference/getGEOSingleCell.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  geoSingleCellManifest("GSE132771")        # GSE: falls back to GSM level
  geoSingleCellManifest("GSM3891612")       # a single sample
} # }
```

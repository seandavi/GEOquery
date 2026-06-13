# Inventory the single-cell supplementary files of a GEO Series

Lists the supplementary files attached to a GSE and classifies each by
single-cell format (10x Matrix Market triplet, 10x HDF5, AnnData h5ad,
loom, Seurat rds, tar archive, or other), extracting the GSM sample id
where present. This lets you see what a single-cell study contains – and
how 10x triplets group by sample – before downloading potentially many
gigabytes.

## Usage

``` r
geoSingleCellManifest(GEO)
```

## Arguments

- GEO:

  A GEO Series accession, e.g. "GSE161228".

## Value

A data.frame with columns `fname`, `sample` (GSM id or NA), `format`,
`role`, and `url`. Zero rows if the GSE has no supplementary files.

## Details

No files are downloaded. The result feeds the planned single-cell
readers (see ADR-0004); reading itself uses Bioconductor importers
(TENxIO, anndataR) that are optional dependencies.

## See also

[`getGEOSuppFiles`](http://seandavi.github.io/GEOquery/reference/getGEOSuppFiles.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  m <- geoSingleCellManifest("GSE161228")
  m
} # }
```

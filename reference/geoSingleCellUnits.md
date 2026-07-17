# Group a single-cell manifest into loadable units

Collapses a
[`geoSingleCellManifest`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md)
into one row per unit and reports completeness. A 10x Matrix Market unit
groups a sample's matrix, barcodes, and features files and is "complete"
only when all three are present; every other format is one unit per
file. The `loadable` column flags units a built-in reader can consume –
complete 10x Matrix Market, 10x HDF5, and AnnData h5ad. loom and Seurat
`.rds` are reported but not loadable (read them with their native
packages).

## Usage

``` r
geoSingleCellUnits(manifest)
```

## Arguments

- manifest:

  A data.frame returned by
  [`geoSingleCellManifest()`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md).

## Value

A data.frame with columns `unit` (the grouping key), `sample`,
`platform` (GPL, or NA), `format`, `n_files`, `status`, and `loadable`.

## See also

[`geoSingleCellManifest`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  m <- geoSingleCellManifest("GSE161228")
  geoSingleCellUnits(m)
} # }
```

# Group a single-cell manifest into loadable units

Collapses a
[`geoSingleCellManifest`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md)
into one row per loadable unit (a sample + format combination) and
reports completeness. A 10x Matrix Market unit is "complete" only when
its matrix, barcodes, and features files are all present; single-file
formats (h5ad, 10x h5, loom, rds) are always complete. The `loadable`
column flags units a reader can consume.

## Usage

``` r
geoSingleCellUnits(manifest)
```

## Arguments

- manifest:

  A data.frame returned by
  [`geoSingleCellManifest()`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md).

## Value

A data.frame with columns `sample`, `format`, `n_files`, `status`, and
`loadable`.

## See also

[`geoSingleCellManifest`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  m <- geoSingleCellManifest("GSE161228")
  geoSingleCellUnits(m)
} # }
```

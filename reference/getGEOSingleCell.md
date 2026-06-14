# Download and read the single-cell data of a GEO Series

High-level, best-effort convenience wrapper: inventories the GSE
([`geoSingleCellManifest`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md)),
groups files into loadable units
([`geoSingleCellUnits`](http://seandavi.github.io/GEOquery/reference/geoSingleCellUnits.md)),
downloads each loadable unit, reads it with
[`readGEOSingleCell`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md),
and returns the results. It reports which units it loads and which it
skips.

## Usage

``` r
getGEOSingleCell(
  GEO,
  samples = NULL,
  format = NULL,
  combine = FALSE,
  destdir = tempdir()
)
```

## Arguments

- GEO:

  A GEO Series accession, e.g. "GSE161228".

- samples:

  Optional character vector of GSM ids to restrict to.

- format:

  Optional format(s) to restrict to ("10x_mtx", "10x_h5", "h5ad").

- combine:

  Logical; if TRUE attempt to `cbind` the per-sample objects into one
  (requires matching features). Default FALSE returns a list.

- destdir:

  Download destination directory.

## Value

A named list of `SingleCellExperiment` (one per sample), or a single
combined object if `combine = TRUE`.

## Details

This handles common, well-structured layouts (clean per-sample 10x or
h5ad). It does NOT handle every GSE: loom and Seurat `.rds` formats,
files packaged inside a `_RAW.tar` archive, and idiosyncratic layouts
(e.g. a single combined matrix for many samples) are out of scope – use
the manifest plus
[`readGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md)
directly for those.

## See also

[`geoSingleCellManifest`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md),
[`readGEOSingleCell`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md)

# Read a single-cell file (or 10x triplet) into a SingleCellExperiment

Low-level reader: given already-downloaded local file(s), dispatch on
format to the appropriate Bioconductor importer and return a
`SingleCellExperiment`. Use this for full control; see
[`getGEOSingleCell`](http://seandavi.github.io/GEOquery/reference/getGEOSingleCell.md)
for the high-level convenience wrapper.

## Usage

``` r
readGEOSingleCell(x, format = NULL)
```

## Arguments

- x:

  A path to a single file (`.h5`/`.h5ad`), a directory containing a 10x
  triplet, or a character vector of the triplet files.

- format:

  One of "10x_mtx", "10x_h5", "h5ad". If NULL (default), guessed from
  `x`.

## Value

A `SingleCellExperiment`.

## Details

Supported formats: `"10x_mtx"` (a directory, or the matrix/barcodes/
features files, read via TENxIO), `"10x_h5"` (CellRanger HDF5, TENxIO),
and `"h5ad"` (AnnData, anndataR). loom and Seurat `.rds` are not
supported here – read them with their native packages.

## See also

[`getGEOSingleCell`](http://seandavi.github.io/GEOquery/reference/getGEOSingleCell.md),
[`geoSingleCellManifest`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md)

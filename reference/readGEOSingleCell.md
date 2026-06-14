# Read a single-cell file (or 10x triplet) into a SingleCellExperiment

Low-level reader: given already-downloaded local file(s), dispatch on
format to the appropriate Bioconductor importer and return a
`SingleCellExperiment` (or a `Seurat` object with `as = "Seurat"`). Use
this for full control; see
[`getGEOSingleCell`](http://seandavi.github.io/GEOquery/reference/getGEOSingleCell.md)
for the high-level convenience wrapper.

## Usage

``` r
readGEOSingleCell(x, format = NULL, as = c("SingleCellExperiment", "Seurat"))
```

## Arguments

- x:

  A path to a single file (`.h5`/`.h5ad`/`.rds`), a directory containing
  a 10x triplet, or a character vector of the triplet files.

- format:

  One of "10x_mtx", "10x_h5", "h5ad", "rds". If NULL (default), guessed
  from `x`.

- as:

  Output class, one of "SingleCellExperiment" (default) or "Seurat"
  (coerced via the Seurat package, an optional dependency).

## Value

A `SingleCellExperiment`, or a `Seurat` object if `as = "Seurat"`.

## Details

Supported formats: `"10x_mtx"` (a directory, or the matrix/barcodes/
features files, read via TENxIO), `"10x_h5"` (CellRanger HDF5, TENxIO),
`"h5ad"` (AnnData, anndataR), and `"rds"` (a saved
`SingleCellExperiment` or `Seurat` object; detected by class). loom is
not supported – read it with `LoomExperiment` directly.

## See also

[`getGEOSingleCell`](http://seandavi.github.io/GEOquery/reference/getGEOSingleCell.md),
[`geoSingleCellManifest`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md)

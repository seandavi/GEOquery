# Download and read the single-cell data of a GEO Series or Sample

High-level, best-effort convenience wrapper: inventories the GSE (or
GSM)
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

  A GEO Series (`"GSE..."`) or Sample (`"GSM..."`) accession, e.g.
  "GSE132771" or "GSM3891612".

- samples:

  Optional character vector of GSM ids to restrict to. Ignored when
  `GEO` is itself a GSM.

- format:

  Optional format(s) to restrict to ("10x_mtx", "10x_h5", "h5ad").

- combine:

  Logical; if TRUE, `cbind` the per-sample objects into one, restricting
  to the features (rownames) common to all samples so they align even
  when samples come from different references or platforms. Errors if
  the samples share no common features (e.g. a study mixing organisms).
  Default FALSE returns a named list.

- destdir:

  Download destination directory.

## Value

A named list of `SingleCellExperiment` (one per sample), or a single
combined object if `combine = TRUE`.

## Details

This handles common, well-structured layouts (clean per-sample 10x or
h5ad), including the very common case where the series ships only a
`_RAW.tar` and the per-sample files live in each GSM suppl directory
(the manifest falls back to the GSM level automatically). You may also
pass a single GSM accession to load just that sample. It does NOT handle
every GSE: loom and Seurat `.rds` formats, files available *only* inside
a `_RAW.tar` archive, and idiosyncratic layouts (e.g. a single combined
matrix for many samples) are out of scope – use the manifest plus
[`readGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md)
directly for those.

## See also

[`geoSingleCellManifest`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md),
[`readGEOSingleCell`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  sce <- getGEOSingleCell("GSM3891612")                  # one sample
  all <- getGEOSingleCell("GSE132771")                   # whole series
  two <- getGEOSingleCell("GSE132771",
                          samples = c("GSM3891612", "GSM3891613"))
} # }
```

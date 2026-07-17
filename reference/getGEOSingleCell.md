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
  by = c("sample", "platform", "all"),
  as = c("SingleCellExperiment", "Seurat"),
  addSampleMeta = TRUE,
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

  Optional format(s) to restrict to ("10x_mtx", "10x_h5", "h5ad",
  "rds").

- by:

  One of `"sample"` (default), `"platform"`, or `"all"` – how to group
  the loaded samples into the return value. See Details.

- as:

  Output class, one of "SingleCellExperiment" (default) or "Seurat"
  (coerced at the boundary via the Seurat package, an optional
  dependency).

- addSampleMeta:

  Logical, default `TRUE`. Attach each sample's GEO metadata – the
  parsed `characteristics_*` fields, `title`, and source name – as
  constant per-cell `colData` columns, so a sample's annotation (age,
  sex, genotype, tissue, treatment, ...) travels with its cells
  (including after combining with `by = "platform"`/`"all"`). Added
  columns are prefixed `"sample."` to avoid colliding with the
  importer's own `colData`. The metadata comes from the Series Matrix
  pData (for a GSE) or the GSM SOFT record (for a lone GSM), so it costs
  one additional small download; whole-study files with no GSM get
  nothing. Metadata lookup is best-effort: a failure warns and is
  skipped, never breaking the data load. Set `FALSE` to skip it.

- destdir:

  Download destination directory.

## Value

Depends on `by`: a named list of objects per sample (`"sample"`); a
named list of combined objects per platform (`"platform"`); or a single
combined object (`"all"`). Each object is a `SingleCellExperiment`, or a
`Seurat` object when `as = "Seurat"`.

## Details

This handles common, well-structured layouts (clean per-sample 10x,
h5ad, or a saved object in `.rds`), including the very common case where
the series ships only a `_RAW.tar` and the per-sample files live in each
GSM suppl directory (the manifest falls back to the GSM level
automatically). You may also pass a single GSM accession to load just
that sample. It does NOT handle every GSE: loom files, files available
*only* inside a `_RAW.tar` archive, and idiosyncratic layouts (e.g. a
single combined matrix for many samples) are out of scope – use the
manifest plus
[`readGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md)
directly for those.

**Grouping (`by`).** A GEO Series has two natural layers – it can span
multiple platforms (GPLs), each holding many samples (GSMs) – and the
platform is the feature-compatibility boundary (samples in one platform
share a feature space; across platforms they generally do not). `by`
chooses the return shape, and the shape is fixed by the argument (not
the data):

- `"sample"` (default):

  a named list with one object per sample (or per whole-study file).

- `"platform"`:

  a named list keyed by platform (GPL), each entry the samples of that
  platform combined into one object. The honest answer for a
  multi-platform study; a single-platform study yields a length-1 list.
  Samples with unknown platform are returned individually.

- `"all"`:

  a single object with every sample combined. Errors if the samples
  share no common features (e.g. a study mixing organisms) – use
  `"platform"` for those.

Combining (for `"platform"`/`"all"`) restricts to the features common to
the group and reconciles per-sample feature annotation so binding
succeeds across CellRanger versions; whole-study single-file formats
already hold one object, so grouping is effectively a no-op for them.

## See also

[`geoSingleCellManifest`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md),
[`readGEOSingleCell`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  sce <- getGEOSingleCell("GSM3891612")                   # one sample
  per_sample <- getGEOSingleCell("GSE132771")             # list by GSM
  per_platform <- getGEOSingleCell("GSE132771", by = "platform")
  # -> list(GPL21103 = <mouse SCE>, GPL24676 = <human SCE>)

  # Each sample's GEO characteristics ride along in colData (default):
  sce <- getGEOSingleCell("GSE125708", by = "all")
  SummarizedExperiment::colData(sce)[, grep("^sample\\.", colnames(
    SummarizedExperiment::colData(sce)))]
  # sample.title, sample.age.ch1, sample.Sex.ch1, sample.tissue.ch1, ...
} # }
```

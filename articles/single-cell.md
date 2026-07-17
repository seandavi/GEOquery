# Single-cell data from GEO

Single-cell studies are now a large fraction of GEO submissions, but
they do not fit the classic GEO mental model. This article explains
*why* single-cell data is awkward to retrieve, the file formats you will
meet, and how GEOquery’s single-cell functions turn them into
`SingleCellExperiment` objects ready for the Bioconductor single-cell
ecosystem.

## Why single-cell data is different

For a microarray or bulk RNA-seq study, the processed matrix lives in
the GEO **Series Matrix** file, and `getGEO("GSE...")` hands you a
`SummarizedExperiment` directly. **Single-cell data almost never works
this way.** The Series Matrix for a single-cell `GSE` is usually empty
or contains only sample-level metadata, because a per-cell matrix with
tens of thousands of columns does not belong in GEO’s sample-by-feature
table.

Instead, the actual data lives in **supplementary files** (see
[Understanding GEO data
formats](http://seandavi.github.io/GEOquery/articles/geo-data-formats.md)).
So a single-cell workflow is fundamentally *supplementary-files-first*:
you must **look at what is attached, decide what is loadable, and then
load it** — which is exactly the shape of the GEOquery single-cell API.

## The file formats you will meet

GEO single-cell submissions are heterogeneous. The common formats:

- **10x Matrix Market triplet** — three files per sample: a sparse
  `matrix.mtx[.gz]`, `barcodes.tsv[.gz]` (cell IDs), and
  `features.tsv[.gz]`/`genes.tsv[.gz]` (gene IDs). The three must be
  read *together*; a missing file makes the sample unreadable.
- **10x HDF5** (`.h5`) — the CellRanger filtered/raw feature-barcode
  matrix in a single HDF5 file.
- **AnnData** (`.h5ad`) — the scverse/Python standard; increasingly
  common on recent GEO submissions.
- **loom**, **Seurat `.rds`** — less common, and intentionally *not*
  handled by GEOquery (read them with `LoomExperiment` / `Seurat`
  directly).

A further wrinkle: files are frequently bundled inside a single
`GSE_RAW.tar` archive, and naming conventions vary wildly between
submitters. This is why the first step is always inspection, not
download.

## Step 1 — inspect: the manifest

[`geoSingleCellManifest()`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md)
lists the supplementary files and classifies each by format and role,
extracting the GSM sample id — **without downloading** anything. It lets
you see, often many gigabytes ahead of time, what a study actually
contains.

``` r

library(GEOquery)
m <- geoSingleCellManifest("GSE161228")
m
#> fname                    sample   format    role       url
#> GSM..._matrix.mtx.gz     GSM...   10x_mtx   matrix     ...
#> GSM..._barcodes.tsv.gz   GSM...   10x_mtx   barcodes   ...
#> GSM..._features.tsv.gz   GSM...   10x_mtx   features   ...
```

## Step 2 — decide: loadable units

[`geoSingleCellUnits()`](http://seandavi.github.io/GEOquery/reference/geoSingleCellUnits.md)
collapses that file list into **loadable units** — one per sample and
format — and tells you whether each is complete. A 10x triplet is only
loadable if all three files are present:

``` r

u <- geoSingleCellUnits(m)
u
#> sample  format   n_files  status                          loadable
#> GSM1    10x_mtx  3        complete                        TRUE
#> GSM2    10x_mtx  2        incomplete (missing features)   FALSE
#> GSM3    h5ad     1        complete                        TRUE
```

This is where you make decisions: which samples, which format if a study
offers more than one, and which incomplete units to skip.

## Step 3 — load

For the common, well-structured cases,
[`getGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/getGEOSingleCell.md)
does the whole thing — manifest → units → download → read — and returns
a **named list of `SingleCellExperiment`, one per sample**. It tells you
what it loads and what it skips, so nothing disappears silently:

``` r

sces <- getGEOSingleCell("GSE161228")
#> Loading GSM1 (10x_mtx)...
#> Loading GSM3 (h5ad)...
#> Skipping 1 unit(s): GSM2 [incomplete (missing features)]

length(sces)      # one SingleCellExperiment per sample
sces[[1]]
```

It returns a list rather than a single combined object on purpose:
per-sample matrices often use different references or feature sets, and
silently reconciling them would be misleading. Combine deliberately when
you know the features match.

### Full control

[`getGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/getGEOSingleCell.md)
deliberately handles only common layouts. For anything unusual — bespoke
naming, a single combined matrix for many samples, files inside a
`_RAW.tar` — use the manifest to find what you want, download with
[`getGEOSuppFiles()`](http://seandavi.github.io/GEOquery/reference/getGEOSuppFiles.md),
and read one unit at a time with
[`readGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md):

``` r

readGEOSingleCell(path_to_dir_or_file)            # format auto-detected
readGEOSingleCell(triplet_files, format = "10x_mtx")
```

## The reader dependencies

Reading uses focused Bioconductor importers, kept as optional
dependencies so a basic GEOquery install stays light:

- **10x mtx / 10x h5** →
  [TENxIO](https://bioconductor.org/packages/TENxIO)
- **h5ad** → [anndataR](https://bioconductor.org/packages/anndataR),
  which reads AnnData natively in R with no Python dependency.

GEOquery will prompt you to install the relevant one if it is missing.

## Downstream: the single-cell ecosystem

Once you have a `SingleCellExperiment`, you are in the heart of
Bioconductor’s single-cell stack. Natural next steps:

- Quality control, normalization, and feature selection with
  [scater](https://bioconductor.org/packages/scater) and
  [scran](https://bioconductor.org/packages/scran).
- Dimensionality reduction, clustering, and annotation following the
  [Orchestrating Single-Cell Analysis with
  Bioconductor](https://bioconductor.org/books/release/OSCA/) (OSCA)
  book.
- On-disk / out-of-memory handling of large matrices with
  [HDF5Array](https://bioconductor.org/packages/HDF5Array) and
  [DelayedArray](https://bioconductor.org/packages/DelayedArray).

## What is intentionally out of scope

GEOquery’s single-cell support targets the discovery-and-load problem,
not everything. It does **not** handle loom or Seurat `.rds` files,
files packaged inside `_RAW.tar`, or idiosyncratic combined-matrix
layouts. The manifest plus
[`readGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md)
is the escape hatch for those, and the design notes are in the project’s
`adr/0004-single-cell-architecture.md`.

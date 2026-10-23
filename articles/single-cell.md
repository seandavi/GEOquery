# Single-cell data from GEO

Single-cell studies are now a large fraction of GEO submissions, but
they do not fit the classic GEO mental model. This article explains
*why* single-cell data is awkward to retrieve, the file formats you will
meet, and how GEOquery’s single-cell functions turn them into
`SingleCellExperiment` objects ready for the Bioconductor single-cell
ecosystem.

``` r

library(GEOquery)
```

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
- **Seurat `.rds`** — a saved Seurat object. GEOquery reads these
  (detected by class) and coerces them to a `SingleCellExperiment`; see
  *Working with Seurat* below.
- **loom** — less common, and intentionally *not* handled by GEOquery
  (read it with `LoomExperiment` directly).

A further wrinkle: files are frequently bundled inside a single
`GSE_RAW.tar` archive, and naming conventions vary wildly between
submitters. This is why the first step is always inspection, not
download.

## Step 1 — inspect: the manifest

[`geoSingleCellManifest()`](http://seandavi.github.io/GEOquery/reference/geoSingleCellManifest.md)
lists the supplementary files and classifies each by format and role,
extracting the GSM sample id and platform — **without downloading**
anything. It lets you see, often many gigabytes ahead of time, what a
study actually contains. GSE132771 is a two-platform study (mouse and
human) where each sample contributes a 10x Matrix Market triplet:

``` r

m <- geoSingleCellManifest("GSE132771")
nrow(m)      # every supplementary file across the study
#> [1] 72
# three files per sample: matrix + barcodes + features (a `url` column, omitted
# here for width, gives each file's download location)
head(m[, c("fname", "sample", "format", "role", "platform")])
#>                                   fname     sample  format     role platform
#> 1 GSM3891612_Bleo1_GFPp_barcodes.tsv.gz GSM3891612 10x_mtx barcodes GPL21103
#> 2    GSM3891612_Bleo1_GFPp_genes.tsv.gz GSM3891612 10x_mtx features GPL21103
#> 3   GSM3891612_Bleo1_GFPp_matrix.mtx.gz GSM3891612 10x_mtx   matrix GPL21103
#> 4 GSM3891613_Bleo2_GFPp_barcodes.tsv.gz GSM3891613 10x_mtx barcodes GPL21103
#> 5    GSM3891613_Bleo2_GFPp_genes.tsv.gz GSM3891613 10x_mtx features GPL21103
#> 6   GSM3891613_Bleo2_GFPp_matrix.mtx.gz GSM3891613 10x_mtx   matrix GPL21103
```

## Step 2 — decide: loadable units

[`geoSingleCellUnits()`](http://seandavi.github.io/GEOquery/reference/geoSingleCellUnits.md)
collapses that file list into **loadable units** — one per sample and
format — and reports whether each is complete and which platform it
belongs to. A 10x triplet is only loadable if all three files
(`matrix`/`barcodes`/`features`) are present; the `status` and
`loadable` columns tell you, so incomplete samples never fail
mid-download:

``` r

u <- geoSingleCellUnits(m)
nrow(u)      # one unit per sample
#> [1] 24
head(u)
#>             unit     sample platform  format n_files   status loadable
#> 1 mtx|GSM3891612 GSM3891612 GPL21103 10x_mtx       3 complete     TRUE
#> 2 mtx|GSM3891613 GSM3891613 GPL21103 10x_mtx       3 complete     TRUE
#> 3 mtx|GSM3891614 GSM3891614 GPL21103 10x_mtx       3 complete     TRUE
#> 4 mtx|GSM3891615 GSM3891615 GPL21103 10x_mtx       3 complete     TRUE
#> 5 mtx|GSM3891616 GSM3891616 GPL21103 10x_mtx       3 complete     TRUE
#> 6 mtx|GSM3891617 GSM3891617 GPL21103 10x_mtx       3 complete     TRUE
```

This is where you make decisions: which samples, which platform, and
which incomplete units (`loadable = FALSE`) to skip.

## Step 3 — load

For the common, well-structured cases,
[`getGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/getGEOSingleCell.md)
does the whole thing — manifest → units → download → read — and returns
a **named list of `SingleCellExperiment`, one per sample**. It tells you
what it loads and what it skips, so nothing disappears silently:

``` r

sces <- getGEOSingleCell("GSE132771")
#> Loading GSM3891612 (10x_mtx)...
#> Loading GSM3891613 (10x_mtx)...
#> ... (one per sample)

length(sces)      # one SingleCellExperiment per sample
sces[[1]]
```

Any unit reported `loadable = FALSE` (an incomplete 10x triplet, say) is
skipped with a message rather than failing the whole call, so nothing
disappears silently.

By default it returns a *list* rather than a single combined object:
per-sample matrices often use different references or feature sets, and
silently reconciling them would be misleading. Because GSE132771 spans
two platforms, combining is meaningful only within a platform —
`by = "platform"` combines the mouse and human samples separately and
refuses to merge across them:

``` r

per_platform <- getGEOSingleCell("GSE132771", by = "platform")
```

### Sample metadata travels with the cells

The single-cell matrices GEO stores as supplementary files carry no
phenotype information — the sample’s age, sex, genotype, tissue,
treatment, and so on live in the GEO **sample record**, not in the
matrix. By default
[`getGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/getGEOSingleCell.md)
fetches that per-sample metadata and attaches it to each object’s
`colData`, broadcast across the sample’s cells, so it is there when you
need it (and it survives combining with `by = "platform"` /
`by = "all"`):

``` r

sces <- getGEOSingleCell("GSE132771")
cd <- SummarizedExperiment::colData(sces[[1]])
cd[, grep("^sample\\.", colnames(cd))]
#>              sample.title sample.source_name_ch1 sample.genotype.ch1 ...
#> AAACCTGAG..  Normal lung  Normal human lung      wild type           ...
#> AAACCTGAGT.. Normal lung  Normal human lung      wild type           ...
```

Added columns are prefixed `sample.` so they never clash with the
importer’s own `colData`. The metadata comes from the Series Matrix (one
small extra download); pass `addSampleMeta = FALSE` to skip it.
Whole-study files that carry no GSM get no sample columns (there is no
single sample to attribute them to).

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

## Working with Seurat

GEOquery speaks `SingleCellExperiment` internally, but Seurat is one
coercion away and supported at both edges (Seurat is an optional
dependency).

Get Seurat objects directly:

``` r

seurat_list <- getGEOSingleCell("GSE161228", as = "Seurat")
# or one object:
seu <- readGEOSingleCell(path, as = "Seurat")
```

Read a Seurat object that a submitter saved as a `.rds` supplementary
file (its contents are detected by class and coerced to a
`SingleCellExperiment`):

``` r

sce <- readGEOSingleCell(path_to_rds)        # Seurat .rds -> SingleCellExperiment
```

Or coerce by hand, in either direction:

``` r

seurat_obj <- Seurat::as.Seurat(sce)
sce <- Seurat::as.SingleCellExperiment(seurat_obj)
```

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
not everything. It does **not** handle loom files, files packaged inside
`_RAW.tar`, or idiosyncratic combined-matrix layouts. The manifest plus
[`readGEOSingleCell()`](http://seandavi.github.io/GEOquery/reference/readGEOSingleCell.md)
is the escape hatch for those, and the design notes are in the project’s
`adr/0004-single-cell-architecture.md`.

## Where to go next

- [Understanding GEO data
  formats](http://seandavi.github.io/GEOquery/articles/geo-data-formats.md)
  — why supplementary files exist and what GEO does and doesn’t parse.
- [From GEO to downstream
  analysis](http://seandavi.github.io/GEOquery/articles/downstream-analysis.md)
  — the OSCA stack from a `SingleCellExperiment`.

# Understanding GEO data formats

This article explains *what* the NCBI Gene Expression Omnibus (GEO)
actually stores, *why* it is shaped the way it is, and *how* GEOquery
maps each form onto a Bioconductor object. Understanding the formats is
the fastest way to predict what
[`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
will return and to debug the occasional surprising record.

## The four GEO entity types

GEO is organized around four accession types, and every GEOquery
workflow starts by recognizing which one you have:

| Prefix | Entity | What it is |
|----|----|----|
| `GPL` | Platform | The array or sequencing platform: the list of probes/features and their annotation. |
| `GSM` | Sample | One hybridization / sequencing run: the measurements for a single biological sample, plus its metadata. |
| `GSE` | Series | A study: a set of related `GSM` samples, optionally spanning several platforms. |
| `GDS` | DataSet | A *curated*, normalized collection assembled by GEO staff from a `GSE`. Far fewer exist; new submissions do not get one. |

A `GSE` is the unit you almost always want. The subtlety is that a
single `GSE` can contain samples from more than one platform, which is
why `getGEO("GSE...")` returns a **list** — one element per platform.

## Two file formats, two code paths

GEO exposes the same underlying data in two very different file formats,
and GEOquery has a distinct parser for each. Knowing which one you are
using explains the object you get back.

### SOFT — the complete, verbose record

SOFT (“Simple Omnibus Format in Text”) is the canonical, loss-less GEO
representation. It is a line-oriented text format where metadata lines
begin with `!`, `^`, or `#`, and data tables are delimited by
`!..._table_begin` / `!..._table_end` markers. SOFT contains
*everything* GEO knows about an entity.

GEOquery parses SOFT into its own S4 classes — `GSE`, `GSM`, `GPL`,
`GDS` — that mirror the file structure. You reach the pieces with
accessors:

``` r

gsm <- getGEO("GSM11805")
Meta(gsm)      # the metadata as a named list
Table(gsm)     # the data table as a data.frame
Columns(gsm)   # descriptions of the data-table columns
```

SOFT is the right choice when you need fields that only exist in the
full record. The cost is size and speed: a large `GSE` family SOFT file
can be enormous, and parsing it is correspondingly slow.

### Series Matrix — the fast, analysis-ready path

For a `GSE`, GEO also publishes a **Series Matrix** file: a compact,
tab-delimited table with the sample metadata as a header block
(`!Sample_...` lines) followed by a single expression matrix (rows =
features, columns = samples). This is what most analyses actually need,
and it parses orders of magnitude faster than SOFT.

This is the GEOquery **default** (`GSEMatrix = TRUE`). The result is a
Bioconductor `SummarizedExperiment` (or, opting back, an `ExpressionSet`
— see below):

``` r

gse <- getGEO("GSE2553")   # returns a list, one SummarizedExperiment per platform
se <- gse[[1]]
assay(se)      # the expression matrix
colData(se)    # sample (phenotype) metadata
rowData(se)    # feature annotation, from the platform GPL
```

If you ever need a field that the Series Matrix omits, drop to SOFT with
`getGEO("GSE...", GSEMatrix = FALSE)`.

### Why `getGEO()` returns different classes

This is the most common source of confusion, and it follows directly
from the formats:

- `GSE` + Series Matrix (default) → a **list** of `SummarizedExperiment`
  (or `ExpressionSet` with `returnType = "ExpressionSet"`)
- `GSE` + SOFT (`GSEMatrix = FALSE`) → a single `GSE` S4 object
- `GSM` / `GPL` / `GDS` → the corresponding S4 object

When in doubt, check [`class()`](https://rdrr.io/r/base/class.html) and
remember that the GSE-matrix path always returns a list because a study
may span platforms.

## ExpressionSet vs. SummarizedExperiment

`SummarizedExperiment` is the modern Bioconductor standard and the
substrate that single-cell (`SingleCellExperiment`) and spatial
(`SpatialExperiment`) classes build on; it is now the GEOquery
**default**. `ExpressionSet` (from **Biobase**) is the historical
container, still available on request:

``` r

# default: SummarizedExperiment
se_list <- getGEO("GSE2553")
# opt back into the legacy ExpressionSet:
eset_list <- getGEO("GSE2553", returnType = "ExpressionSet")
# or convert an existing ExpressionSet without re-downloading:
se <- as_SummarizedExperiment(eset_list[[1]])
```

The accessor vocabulary differs —
[`assay()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)/[`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)/[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
for `SummarizedExperiment` versus
[`exprs()`](https://rdrr.io/pkg/Biobase/man/exprs.html)/[`pData()`](https://rdrr.io/pkg/Biobase/man/phenoData.html)/[`fData()`](https://rdrr.io/pkg/Biobase/man/featureData.html)
for `ExpressionSet` — so pick one and stay consistent.

## Supplementary files: where everything else lives

Series Matrix and SOFT cover *processed* expression tables, but a great
deal of GEO data — raw reads, processed counts, single-cell matrices,
peak calls, images — is attached as **supplementary files** that GEO
does not parse and whose format it does not constrain. GEOquery lists
and downloads them but, by design, does not try to interpret arbitrary
content:

``` r

# see what is attached, without downloading
getGEOSuppFiles("GSE63137", fetch_files = FALSE)
```

This matters enormously for sequencing-era data. RNA-seq counts and
**all** single-cell data live here, not in the Series Matrix — which is
why those have dedicated entry points (see the
[RNA-seq](http://seandavi.github.io/GEOquery/articles/rnaseq.md) and
[single-cell](http://seandavi.github.io/GEOquery/articles/single-cell.md)
articles).

## Where to go next

- [Single-cell data from
  GEO](http://seandavi.github.io/GEOquery/articles/single-cell.md) — why
  single-cell lives in supplementary files and how to load it into a
  `SingleCellExperiment`.
- [Downstream
  analysis](http://seandavi.github.io/GEOquery/articles/downstream-analysis.md)
  — taking a GEOquery object into differential expression and beyond,
  with pointers to the relevant Bioconductor packages.

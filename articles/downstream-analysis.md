# From GEO to downstream analysis

GEOquery’s job ends once your data is a Bioconductor object. This
article is a map of *where to go next* — how the object you got from
[`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
connects to the rest of the Bioconductor ecosystem — rather than a
tutorial for any one method. The goal is to save you the “I have a
`SummarizedExperiment`, now what?” moment.

## Know your object

What you do next depends on what GEOquery returned (see [Understanding
GEO data
formats](http://seandavi.github.io/GEOquery/articles/geo-data-formats.md)):

| You have | Typical source | Downstream entry point |
|----|----|----|
| `SummarizedExperiment` (list) | GSE Series Matrix (default), RNA-seq counts | limma / DESeq2 / edgeR |
| `ExpressionSet` (list) | GSE Series Matrix with `returnType = "ExpressionSet"` | limma |
| `SingleCellExperiment` | single-cell readers | scater / scran / OSCA |
| `GSE`/`GSM`/`GPL`/`GDS` S4 | SOFT parsing | accessors, then convert |

## Microarray: limma

Microarray Series Matrix data arrives already processed (typically
log-transformed, normalized intensities). The standard path is a linear
model with [limma](https://bioconductor.org/packages/limma), which
accepts a matrix:

``` r

library(limma)
library(SummarizedExperiment)
se <- getGEO("GSE2553")[[1]]                 # SummarizedExperiment (default)

# `group` here stands in for a real grouping column you derive from the sample
# metadata -- see the note below; it will not exist verbatim in colData(se).
design <- model.matrix(~ group, data = colData(se))
fit <- eBayes(lmFit(assay(se), design))
topTable(fit, coef = 2)
```

The hardest part is usually not the model but extracting clean grouping
variables from `colData(se)` — GEO sample metadata is free text, so
expect to parse `characteristics_ch1` fields into the factor you
actually model.

## RNA-seq counts: DESeq2 / edgeR / limma-voom

For NCBI-computed RNA-seq counts (see the [RNA-seq
article](http://seandavi.github.io/GEOquery/articles/rnaseq.md)), use a
count-based model. The `SummarizedExperiment` GEOquery returns plugs
directly into [DESeq2](https://bioconductor.org/packages/DESeq2):

``` r

library(DESeq2)
se <- getRNASeqData("GSE164073")
# `condition` must be a real column you added to colData(se) first.
dds <- DESeqDataSet(se, design = ~ condition)
dds <- DESeq(dds)
results(dds)
```

[edgeR](https://bioconductor.org/packages/edgeR) and
[limma-voom](https://bioconductor.org/packages/limma) are equally good
choices and consume the same count matrix.

## Single cell: the OSCA stack

A `SingleCellExperiment` from the [single-cell
readers](http://seandavi.github.io/GEOquery/articles/single-cell.md) is
the entry point to Bioconductor’s single-cell ecosystem —
[scater](https://bioconductor.org/packages/scater),
[scran](https://bioconductor.org/packages/scran), and the [OSCA
book](https://bioconductor.org/books/release/OSCA/) for the full
quality-control → normalization → clustering → annotation arc.

## Converting and annotating

Two recurring needs:

- **Modernize a legacy result.**
  [`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
  returns `SummarizedExperiment` by default now, but if you have an
  older `ExpressionSet` (or asked for one with
  `returnType = "ExpressionSet"`), convert it without re-downloading:

  ``` r

  se <- as_SummarizedExperiment(getGEO("GSE2553", returnType = "ExpressionSet")[[1]])
  ```

- **Re-annotate features.** GEO platform annotation can be dated. For
  mapping probe/gene identifiers, reach for the Bioconductor annotation
  infrastructure —
  [AnnotationDbi](https://bioconductor.org/packages/AnnotationDbi),
  organism packages such as `org.Hs.eg.db`, and
  [biomaRt](https://bioconductor.org/packages/biomaRt).

## Reproducibility

GEO records can change and large downloads are slow, so cache
deliberately and record versions. The [Finding and downloading
data](http://seandavi.github.io/GEOquery/articles/finding-and-downloading-data.md)
article covers persistent caching (`destdir=` and the `BiocFileCache`
layer); alongside that, capture
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) with your
results so the GEOquery and dependency versions are pinned.

## Where to go next

- [Finding and downloading
  data](http://seandavi.github.io/GEOquery/articles/finding-and-downloading-data.md)
  — caching and reproducible retrieval.
- [RNA-seq](http://seandavi.github.io/GEOquery/articles/rnaseq.md) and
  [Single-cell](http://seandavi.github.io/GEOquery/articles/single-cell.md)
  — the data-type paths that feed these workflows.

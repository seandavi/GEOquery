# RNA-seq quantifications from GEO

A major barrier to reusing public RNA-seq data has been that every
submitter processes raw reads differently — different aligners,
references, and counting rules — so counts from two studies are rarely
directly comparable. To address this, NCBI computes a **uniform set of
RNA-seq quantifications** for many GEO Series using a single consistent
pipeline. This article explains what those are and how GEOquery
retrieves them.

``` r

library(GEOquery)
```

## NCBI-computed quantifications vs. submitter files

There are two distinct things you might mean by “the counts” for a GEO
RNA-seq study:

1.  **Submitter-provided processed files** — whatever the authors
    uploaded as supplementary files. Formats and gene models vary by
    study; GEOquery lists and downloads them
    ([`getGEOSuppFiles()`](http://seandavi.github.io/GEOquery/reference/getGEOSuppFiles.md))
    but does not interpret them.
2.  **NCBI-computed quantifications** — a raw-count matrix (and
    accompanying annotation) that NCBI generated uniformly from the
    study’s SRA runs. These are comparable *across* studies and are what
    the GEOquery RNA-seq functions target.

Not every Series has NCBI quantifications; you can check first.

## Checking availability and loading

``` r

# Does this Series have NCBI-computed quantifications?
hasRNASeqQuantifications("GSE164073")
#> [1] TRUE
```

``` r

# Load them as a SummarizedExperiment
se <- getRNASeqData("GSE164073")
se
#> class: SummarizedExperiment 
#> dim: 39376 18 
#> metadata(5): experimentData annotation protocolData genomeInfo
#>   created_at
#> assays(1): counts
#> rownames(39376): 100287102 653635 ... 4576 4571
#> rowData names(18): GeneID Symbol ... GOProcess GOComponent
#> colnames(18): GSM4996084 GSM4996085 ... GSM4996100 GSM4996101
#> colData names(45): title geo_accession ... time.point.ch1 tissue.ch1
```

The result is a `SummarizedExperiment`: `assay(se)` is the raw count
matrix, `rowData(se)` carries gene annotation, and `colData(se)` the
sample metadata. Because it is a `SummarizedExperiment`, it drops
straight into the differential expression workflows described in the
[downstream
analysis](http://seandavi.github.io/GEOquery/articles/downstream-analysis.md)
article.

``` r

assay(se)[1:5, 1:3]
#>           GSM4996084 GSM4996085 GSM4996086
#> 100287102          2          5          3
#> 653635           244        236        337
#> 102466751         25         17         34
#> 107985730          1          1          1
#> 100302278          0          0          0
```

## Why raw counts (and what to do with them)

The NCBI quantifications are deliberately **raw counts**, not normalized
values, because the correct normalization depends on your downstream
method. Count-based differential expression tools expect raw counts and
model the counting noise themselves:

- [DESeq2](https://bioconductor.org/packages/DESeq2)
- [edgeR](https://bioconductor.org/packages/edgeR)
- [limma](https://bioconductor.org/packages/limma) with `voom`

See [downstream
analysis](http://seandavi.github.io/GEOquery/articles/downstream-analysis.md)
for a worked sketch.

## When you need the submitter’s files instead

If a study has no NCBI quantifications, or you specifically need the
authors’ processed output (e.g. TPM, or a custom gene model), fall back
to the supplementary files and read them yourself:

``` r

getGEOSuppFiles("GSE164073", fetch_files = FALSE)   # inspect first
#>                               fname
#> 1 GSE164073_Eye_count_matrix.csv.gz
#>                                                                                                   url
#> 1 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164073/suppl/GSE164073_Eye_count_matrix.csv.gz
```

## Where to go next

- [Finding and downloading
  data](http://seandavi.github.io/GEOquery/articles/finding-and-downloading-data.md)
  — search for RNA-seq studies and cache large downloads.
- [From GEO to downstream
  analysis](http://seandavi.github.io/GEOquery/articles/downstream-analysis.md)
  — DESeq2 / edgeR / limma-voom from a count `SummarizedExperiment`.

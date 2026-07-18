# Getting started with GEOquery

GEOquery is the bridge between the NCBI Gene Expression Omnibus (GEO)
and Bioconductor: it downloads a GEO accession and parses it into a
ready-to-use Bioconductor object. This page gets you from install to a
first analysis-ready object; the other vignettes go deeper on formats,
finding data, and specific data types (links at the end).

## Install

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("GEOquery")
```

## Your first record

Everything starts with
[`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md) and
a GEO accession. The most common case is a GEO **Series** (a `GSE` — one
study), which
[`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
returns as a list of `SummarizedExperiment` objects, one per platform:

``` r

library(GEOquery)

gse <- getGEO("GSE2553")
length(gse)          # one element per platform in the study
#> [1] 1
se <- gse[[1]]
se
#> class: RangedSummarizedExperiment 
#> dim: 12600 181 
#> metadata(3): experimentData annotation protocolData
#> assays(1): exprs
#> rownames(12600): 1 2 ... 12599 12600
#> rowData names(13): ID PenAt ... LLID Chimeric_Cluster_IDs
#> colnames(181): GSM48681 GSM48682 ... GSM48860 GSM48861
#> colData names(30): title geo_accession ... supplementary_file
#>   data_row_count
```

A study is returned as a *list* because a single GEO Series can span
more than one platform. Here there is just one, so we take `gse[[1]]`.

The `SummarizedExperiment` bundles three aligned pieces — the expression
matrix, the per-sample metadata, and the per-feature annotation —
reached with three accessors:

``` r

assay(se)[1:5, 1:3]                       # expression matrix (features x samples)
#>     GSM48681   GSM48682   GSM48683
#> 1  0.2701103  0.3925373  0.4186763
#> 2  6.3459203  2.1304703  2.0750533
#> 3 -0.0918793 -0.2411003 -0.2499943
#> 4  1.4679053  0.6252703  0.4673843
#> 5 -0.2817373  0.6492023 -1.4973643
```

``` r

colData(se)[1:3, c("title", "geo_accession", "source_name_ch1")]   # sample metadata
#> DataFrame with 3 rows and 3 columns
#>                           title geo_accession     source_name_ch1
#>                     <character>   <character>         <character>
#> GSM48681 Patient sample ST18,..      GSM48681 Dermatofibrosarcoma
#> GSM48682 Patient sample ST410..      GSM48682       Ewing Sarcoma
#> GSM48683 Patient sample ST130..      GSM48683        Sarcoma, NOS
```

``` r

head(rowData(se))                         # feature annotation, from the platform
#> DataFrame with 6 rows and 13 columns
#>          ID     PenAt     RowAt  ColumnAt      CLONE_ID     SPOT_ID    PlatePos
#>   <integer> <integer> <integer> <integer>   <character> <character> <character>
#> 1         1         1         1         1  IMAGE:502055                 HsKG1A1
#> 2         2         1         1         2  IMAGE:511814                HsKG20A1
#> 3         3         1         1         3   IMAGE:79592                HsKG40A1
#> 4         4         1         1         4 IMAGE:1571993                HsKG60A1
#> 5         5         1         1         5  IMAGE:150221                HsKG84A1
#> 6         6         1         1         6  IMAGE:591632                FHsKG9A1
#>       UNIGENE                   Name      Symbol                Aliases
#>   <character>            <character> <character>            <character>
#> 1   Hs.149103        Arylsulfatase B        ARSB ||MPS6||ARSB||ASB||G..
#> 2   Hs.435302 Zinc finger protein ..        ZNF3 ||ZNF3||zinc finger ..
#> 3   Hs.512807 Aldo-keto reductase ..      AKR7A2 ||AKR7A2||AFAR||AFLA..
#> 4    Hs.11590            Cathepsin F        CTSF  ||CTSF||cathepsin F||
#> 5   Hs.310645 RAB1A, member RAS on..       RAB1A ||RAB1A||RAS-ASSOCIA..
#> 6   Hs.460317 Amyotrophic lateral ..        ALS4 ||ALS4||Amyotrophic ..
#>          LLID Chimeric_Cluster_IDs
#>   <character>          <character>
#> 1         411         Not chimeric
#> 2        7551         Not chimeric
#> 3        8574         Not chimeric
#> 4        8722         Not chimeric
#> 5        5861         Not chimeric
#> 6       23064         Not chimeric
```

That is the whole loop: one accession in, an analysis-ready object out.
(Prefer the legacy `ExpressionSet`? Pass
`returnType = "ExpressionSet"`.)

## The other entity types

GEO has four accession types. A `GSE` (Series) is the one you usually
want, but Samples (`GSM`), Platforms (`GPL`), and curated DataSets
(`GDS`) parse to GEOquery’s own S4 classes:

``` r

class(getGEO("GSM11805"))     # a single sample
#> [1] "GSM"
#> attr(,"package")
#> [1] "GEOquery"
class(getGEO("GDS507"))       # a curated dataset
#> [1] "GDS"
#> attr(,"package")
#> [1] "GEOquery"
```

What these classes are, and why
[`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
returns different things for different accessions, is the subject of
[Understanding GEO data
formats](http://seandavi.github.io/GEOquery/articles/geo-data-formats.md).

## Peeking at supplementary files

Processed expression tables are only part of a study. Raw data, RNA-seq
counts, and single-cell matrices arrive as **supplementary files**. You
can list them without downloading anything:

``` r

getGEOSuppFiles("GSE63137", fetch_files = FALSE)
#>                                                               fname
#> 1                   GSE63137_ATAC-seq_PV_neurons_HOMER_peaks.bed.gz
#> 2                  GSE63137_ATAC-seq_VIP_neurons_HOMER_peaks.bed.gz
#> 3           GSE63137_ATAC-seq_excitatory_neurons_HOMER_peaks.bed.gz
#> 4   GSE63137_ChIP-seq_H3K27ac_excitatory_neurons_SICER_peaks.bed.gz
#> 5  GSE63137_ChIP-seq_H3K27me3_excitatory_neurons_SICER_peaks.bed.gz
#> 6   GSE63137_ChIP-seq_H3K4me1_excitatory_neurons_SICER_peaks.bed.gz
#> 7   GSE63137_ChIP-seq_H3K4me3_excitatory_neurons_SICER_peaks.bed.gz
#> 8                         GSE63137_MethylC-seq_DMRs_methylpy.txt.gz
#> 9                  GSE63137_MethylC-seq_PV_neurons_UMRs_LMRs.txt.gz
#> 10                GSE63137_MethylC-seq_VIP_neurons_UMRs_LMRs.txt.gz
#> 11         GSE63137_MethylC-seq_excitatory_neurons_UMRs_LMRs.txt.gz
#> 12                                                 GSE63137_RAW.tar
#>                                                                                                                                 url
#> 1                   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_ATAC-seq_PV_neurons_HOMER_peaks.bed.gz
#> 2                  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_ATAC-seq_VIP_neurons_HOMER_peaks.bed.gz
#> 3           https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_ATAC-seq_excitatory_neurons_HOMER_peaks.bed.gz
#> 4   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_ChIP-seq_H3K27ac_excitatory_neurons_SICER_peaks.bed.gz
#> 5  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_ChIP-seq_H3K27me3_excitatory_neurons_SICER_peaks.bed.gz
#> 6   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_ChIP-seq_H3K4me1_excitatory_neurons_SICER_peaks.bed.gz
#> 7   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_ChIP-seq_H3K4me3_excitatory_neurons_SICER_peaks.bed.gz
#> 8                         https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_MethylC-seq_DMRs_methylpy.txt.gz
#> 9                  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_MethylC-seq_PV_neurons_UMRs_LMRs.txt.gz
#> 10                https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_MethylC-seq_VIP_neurons_UMRs_LMRs.txt.gz
#> 11         https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_MethylC-seq_excitatory_neurons_UMRs_LMRs.txt.gz
#> 12                                                 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl/GSE63137_RAW.tar
```

## Where to go next

- [Understanding GEO data
  formats](http://seandavi.github.io/GEOquery/articles/geo-data-formats.md)
  — the four entity types, SOFT vs. Series Matrix, and why
  [`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
  returns different classes.
- [Finding and downloading
  data](http://seandavi.github.io/GEOquery/articles/finding-and-downloading-data.md)
  — search GEO from R, control what
  [`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
  fetches, cache downloads, and reach private records.
- [RNA-seq quantifications from
  GEO](http://seandavi.github.io/GEOquery/articles/rnaseq.md) — NCBI’s
  uniformly-computed counts.
- [Single-cell data from
  GEO](http://seandavi.github.io/GEOquery/articles/single-cell.md) — the
  inspect → decide → load workflow into a `SingleCellExperiment`.
- [From GEO to downstream
  analysis](http://seandavi.github.io/GEOquery/articles/downstream-analysis.md)
  — taking a GEOquery object into limma / DESeq2 / edgeR / the
  single-cell ecosystem.

## Getting help

- Usage questions: the [Bioconductor support
  site](https://support.bioconductor.org/), tagged `geoquery`.
- Bugs and feature requests: the [issue
  tracker](https://github.com/seandavi/GEOquery/issues) — please include
  a GEO accession and
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).

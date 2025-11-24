# Get GEO RNA-seq quantifications as a SummarizedExperiment object

For human and mouse GEO datasets, NCBI GEO attempts to process the raw
data and provide quantifications in the form of raw counts and an
annotation file. This function downloads the raw counts and annotation
files from GEO and merges that with the metadata from the GEO object to
create a SummarizedExperiment.

## Usage

``` r
getRNASeqData(accession)
```

## Arguments

- accession:

  GEO accession number

## Value

A SummarizedExperiment object with the raw counts as the counts assay,
the annotation as the rowData, and the metadata from GEO as the colData.

## Details

A major barrier to fully exploiting and reanalyzing the massive volumes
of public RNA-seq data archived by SRA is the cost and effort required
to consistently process raw RNA-seq reads into concise formats that
summarize the expression results. To help address this need, the NCBI
SRA and GEO teams have built a pipeline that precomputes RNA-seq gene
expression counts and delivers them as count matrices that may be
incorporated into commonly used differential expression analysis and
visualization software.

The pipeline processes RNA-seq data from SRA using the HISAT2 aligner
and and then generates gene expression counts using the featureCounts
program.

See the [GEO
documentation](https://ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html) for
more details.

## Examples

``` r
se <- getRNASeqData("GSE164073")
#> Found 1 file(s)
#> GSE164073_series_matrix.txt.gz
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

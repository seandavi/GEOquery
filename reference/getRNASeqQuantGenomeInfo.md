# Extract genome build and species for GEO RNA-seq quantification

This function extracts the genome build and species information for a
GEO RNA-seq quantification.

## Usage

``` r
getRNASeqQuantGenomeInfo(gse)
```

## Arguments

- gse:

  GEO accession number

## Value

A character vector with the genome build and species information

## Examples

``` r
getRNASeqQuantGenomeInfo("GSE164073")
#>                    genome_build                         species 
#>                    "GRCh38.p13"                         "Human" 
#>                           fname 
#> "Human.GRCh38.p13.annot.tsv.gz" 
```

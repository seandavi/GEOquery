# Get RNA-seq quantification and annotation from GEO

This function downloads the raw counts and annotation files from GEO for
a given GEO accession number.

## Usage

``` r
getRNASeqQuantResults(gse)
```

## Arguments

- gse:

  GEO accession number

## Value

A list with two elements: quants (a matrix of raw counts) and annotation
(a data frame of annotation information).

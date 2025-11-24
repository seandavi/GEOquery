# Read raw counts from GEO

This function reads the raw counts from a GEO link. The raw counts are
expected to be in a tab-separated file with the first column containing
the gene IDs and the remaining columns containing the raw counts.

## Usage

``` r
readRNAQuantRawCounts(link)
```

## Arguments

- link:

  A link to the raw counts file

## Value

A matrix of raw counts with gene IDs as row names

## Details

This function reads the raw counts and returns a matrix with the gene
IDs as the row names, ready for use in creating a SummarizedExperiment.

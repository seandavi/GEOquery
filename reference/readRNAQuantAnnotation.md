# Read RNA-seq quantification annotation from GEO

This function reads the annotation file from a GEO link. The annotation
file is expected to be a tab-separated file with the first column
containing the gene IDs and the remaining columns containing the
annotation information.

## Usage

``` r
readRNAQuantAnnotation(link)
```

## Arguments

- link:

  A link to the annotation file

## Value

A data frame of annotation information with gene IDs as row names

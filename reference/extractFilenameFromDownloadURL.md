# Extract filename from a GEO download URL

This function extracts the filename from a GEO download URL. The
filename is expected to be a query parameter called "file". If the query
parameter is not found, the function returns NULL.

## Usage

``` r
extractFilenameFromDownloadURL(url)
```

## Arguments

- url:

  A GEO download URL

## Value

A character vector with the filename

## Details

The idea is to use this function to extract filenames that contain
important metadata from the GEO RNA-seq quantification.

In particular, the filename is expected to contain the genome build and
species information that we can attach to the SummarizedExperiment.

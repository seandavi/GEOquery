# GSE Supplemental file listing

The GEO Series records often have one or more supplemental files. In
most cases, those files are archived as '.tar' files, the contents of
which are only available in a file listing file not present on the
website for download.

## Usage

``` r
getGEOSeriesFileListing(GSE)
```

## Arguments

- GSE:

  character(1) the GSE accession

## Value

A data.frame with 5 columns. See example.

## Details

This function reads that file listing file and returns the results as a
data.frame.

## Examples

``` r
if (FALSE) { # \dontrun{
getGEOSeriesFileListing('GSE288770')

} # }
```

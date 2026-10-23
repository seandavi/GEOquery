# Get GEO supplemental file URL for a given GEO accession

Get GEO supplemental file URL for a given GEO accession

## Usage

``` r
getGEOSuppFileURL(GEO)
```

## Arguments

- GEO:

  A GEO accession, e.g. a Series ("GSE..."), Sample ("GSM..."), or
  Platform ("GPL...") id. The accession type determines which FTP
  `suppl/` directory URL is returned.

## Value

A character(1) URL of the accession's supplementary-file directory on
the NCBI GEO FTP site.

## Examples

``` r
# an example of a GEO supplemental file URL
# with a set of single-cell RNA-seq data
url = getGEOSuppFileURL("GSE161228")
url
#> [1] "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161228/suppl/"

if (FALSE) { # \dontrun{
  browseURL(url)
} # }
```

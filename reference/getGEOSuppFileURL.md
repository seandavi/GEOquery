# Get GEO supplemental file URL for a given GEO accession

Get GEO supplemental file URL for a given GEO accession

## Usage

``` r
getGEOSuppFileURL(GEO)
```

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

# The URL for a GEO accession

Sometimes, you just need the URL for a GEO accession. This function
returns the URL for a given GEO accession number that can be used to
access the GEO page for that accession.

## Usage

``` r
urlForAccession(geo)
```

## Arguments

- geo:

  A GEO accession number

## Value

A character vector with the URL for the GEO accession

## See also

[`browseGEOAccession`](http://seandavi.github.io/GEOquery/reference/browseGEOAccession.md)

## Examples

``` r
urlForAccession("GSE262484")
#> [1] "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262484"
```

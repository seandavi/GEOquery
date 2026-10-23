# Does a GEO accession have RNA-seq quantifications?

This function checks if a GEO accession number has RNA-seq
quantifications available. It does this by checking if the GEO accession
number has a "RNA-Seq raw counts" link available on the GEO download
page.

## Usage

``` r
hasRNASeqQuantifications(accession)
```

## Arguments

- accession:

  GEO accession number

## Value

TRUE if the GEO accession number has RNA-seq quantifications available,
FALSE otherwise.

## Examples

``` r
hasRNASeqQuantifications("GSE164073")
#> [1] TRUE
```

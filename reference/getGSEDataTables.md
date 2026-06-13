# Get GSE data tables from GEO into R data structures.

In some cases, instead of individual sample records (GSM) containing
information regarding sample phenotypes, the GEO Series contains that
information in an attached data table. And example is given by GSE3494
where there are two data tables with important information contained
within them. Using getGEO with the standard parameters downloads the
GSEMatrix file which, unfortunately, does not contain the information in
the data tables. This function simply downloads the “header” information
from the GSE record and parses out the data tables into R data.frames.

## Usage

``` r
getGSEDataTables(GSE)
```

## Arguments

- GSE:

  The GSE identifier, such as “GSE3494”.

## Value

A list of data.frames.

## See also

[`getGEO`](http://seandavi.github.io/GEOquery/reference/getGEO.md)

## Author

Sean Davis <sdavis2@mail.nih.gov>

## Examples

``` r
if (FALSE) { # \dontrun{

dfl = getGSEDataTables('GSE3494')
lapply(dfl,head)


} # }
```

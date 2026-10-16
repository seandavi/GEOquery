# Get GSE data tables from GEO into R data structures.

In some cases, instead of individual sample records (GSM) containing
information regarding sample phenotypes, the GEO Series contains that
information in one or more data tables attached at the Series level. An
example is given by GSE3494, where there are two data tables with
important information contained within them. Series-level per-cell
annotation tables from single-cell studies (the `!series_table` blocks
in SOFT format, e.g. the "Listing of Individual Cells" table in
GSE98638) are exposed here as well. Using getGEO with the standard
parameters downloads the GSEMatrix file which, unfortunately, does not
contain the information in the data tables. This function simply
downloads the “header” information from the GSE record and parses out
the data tables into R data.frames.

## Usage

``` r
getGSEDataTables(GSE)
```

## Arguments

- GSE:

  The GSE identifier, such as “GSE3494”.

## Value

A list of data.frames, one per `<Data-Table>` in the Series record (a
Series may carry zero, one, or several). Each data.frame's column names
are taken from the table's column definitions.

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

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

dfl = getGSEDataTables('GSE3494')
#> Warning: Failed to open 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&form=xml&view=full&acc=GSE3494': The requested URL returned error: 504
#> Error in open.connection(x, "rb"): cannot open the connection
lapply(dfl,head)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'X' in selecting a method for function 'lapply': object 'dfl' not found

```

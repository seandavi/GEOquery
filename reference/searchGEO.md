# Search GEO database

This function searches the [GDS](https://www.ncbi.nlm.nih.gov/gds)
database, and return a data.frame for all the search results.

## Usage

``` r
searchGEO(query, step = 500L)
```

## Arguments

- query:

  character, the search term. The NCBI uses a search term syntax which
  can be associated with a specific search field with square brackets.
  So, for instance "Homo sapiens\[ORGN\]" denotes a search for
  `Homo sapiens` in the “Organism” field. Details see
  <https://www.ncbi.nlm.nih.gov/geo/info/qqtutorial.html>. The names and
  definitions of these fields can be identified using
  [searchFieldsGEO](http://seandavi.github.io/GEOquery/reference/searchFieldsGEO.md).

- step:

  the number of records to fetch from the database each time. You may
  choose a smaller value if failed.

## Value

a data.frame contains the search results

## Details

The NCBI allows users to access more records (10 per second) if they
register for and use an API key.
[set_entrez_key](https://rdrr.io/pkg/rentrez/man/set_entrez_key.html)
function allows users to set this key for all calls to rentrez functions
during a particular R session. You can also set an environment variable
`ENTREZ_KEY` by [Sys.setenv](https://rdrr.io/r/base/Sys.setenv.html).
Once this value is set to your key rentrez will use it for all requests
to the NCBI. Details see
<https://docs.ropensci.org/rentrez/articles/rentrez_tutorial.html#rate-limiting-and-api-keys>

## See also

[searchFieldsGEO](http://seandavi.github.io/GEOquery/reference/searchFieldsGEO.md)

## Examples

``` r
if (FALSE) { # \dontrun{
searchGEO("diabetes[ALL] AND Homo sapiens[ORGN] AND GSE[ETYP]")
} # }
```

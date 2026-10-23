# Download a file from GEO soft file to the local machine

This function simply downloads a SOFT format file associated with the
GEO accession number given.

## Usage

``` r
getGEOfile(
  GEO,
  destdir = tempdir(),
  AnnotGPL = FALSE,
  amount = c("full", "brief", "quick", "data")
)
```

## Arguments

- GEO:

  Character string, the GEO accession for download (eg., GDS84, GPL96,
  GSE2553, or GSM10)

- destdir:

  Directory in which to store the resulting downloaded file. Defaults to
  tempdir()

- AnnotGPL:

  A boolean defaulting to FALSE as to whether or not to use the
  Annotation GPL information. These files are nice to use because they
  contain up-to-date information remapped from Entrez Gene on a regular
  basis. However, they do not exist for all GPLs; in general, they are
  only available for GPLs referenced by a GDS

- amount:

  Amount of information to pull from GEO. Only applies to GSE, GPL, or
  GSM. See details...

## Value

Invisibly returns the full path of the downloaded file.

## Details

This function downloads GEO SOFT files based on accession number. It
does not do any parsing. The first two arguments should be fairly
self-explanatory, but the last is based on the input to the acc.cgi url
at the geo website. In the default 'full' mode, the entire SOFT format
file is downloaded. Both 'brief' and 'quick' offer shortened versions of
the files, good for 'peeking' at the file before a big download on a
slow connection. Finally, 'data' downloads only the data table part of
the SOFT file and is good for downloading a simple EXCEL-like file for
use with other programs (a convenience).

## References

http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

## See also

[`getGEO`](http://seandavi.github.io/GEOquery/reference/getGEO.md)

## Author

Sean Davis

## Examples

``` r

# myfile <- getGEOfile('GDS10')
```

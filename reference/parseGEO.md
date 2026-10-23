# Parse GEO text

Workhorse GEO parsers.

## Usage

``` r
parseGEO(
  fname,
  GSElimits,
  destdir = tempdir(),
  AnnotGPL = FALSE,
  getGPL = TRUE,
  parseCharacteristics = TRUE
)
```

## Arguments

- fname:

  The filename of a SOFT format file. If the filename ends in .gz, a
  gzfile() connection is used to read the file directly.

- GSElimits:

  Used to limit the number of GSMs parsed into the GSE object; useful
  for memory management for large GSEs.

- destdir:

  The destination directory into which files will be saved (to be used
  for caching)

- AnnotGPL:

  Fetch the annotation GPL if available

- getGPL:

  Fetch the GPL associated with a GSEMatrix entity (should remain TRUE
  for all normal use cases)

- parseCharacteristics:

  Whether or not to parse the characteristics information (if available)
  for a GSE Matrix file. Set to FALSE if you experience trouble parsing
  the characteristics.

## Value

parseGEO returns an object of the associated type. For example, if it is
passed the text from a GDS entry, a GDS object is returned.

## Details

These are probably not useful to the end-user. Use getGEO to access
these functions. parseGEO simply delegates to the appropriate specific
parser. There should be no reason to use the parseGPL, parseGDS,
parseGSE, or parseGSM functions directly.

## See also

[`getGEO`](http://seandavi.github.io/GEOquery/reference/getGEO.md)

## Author

Sean Davis

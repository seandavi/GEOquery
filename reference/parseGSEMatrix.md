# Parse a GSE mstrix file

Not meant for user calling, but parses a single GSEMatrix file.

## Usage

``` r
parseGSEMatrix(
  fname,
  AnnotGPL = FALSE,
  destdir = tempdir(),
  getGPL = TRUE,
  parseCharacteristics = TRUE
)
```

## Arguments

- fname:

  filename

- AnnotGPL:

  set to TRUE to get the annotation GPL version

- destdir:

  the destdination directory for download

- getGPL:

  whether or not to get the GPL associated

- parseCharacteristics:

  Whether or not to do full 'characteristic' parsing

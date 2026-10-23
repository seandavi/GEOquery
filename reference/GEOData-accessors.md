# Accessors for GEOquery objects

Accessor generics for the S4 objects returned by
[`getGEO`](http://seandavi.github.io/GEOquery/reference/getGEO.md) when
parsing SOFT-format records (`GSE`, `GSM`, `GPL`, `GDS`). Use these
rather than reaching into slots directly.

## Arguments

- object:

  A GEOquery S4 object (`GSE`, `GSM`, `GPL`, `GDS`, or `GEODataTable`).

## Value

`Meta()` a list; `Accession()` a character string; `Table()` and
`Columns()` data.frames; `dataTable()` a `GEODataTable`; `GSMList()` and
`GPLList()` named lists.

## Details

- `Meta(object)`:

  The record metadata as a named list (title, submission dates,
  sample/platform attributes, and so on).

- `Accession(object)`:

  The GEO accession (the `geo_accession` metadata field).

- `Table(object)`:

  The data table as a `data.frame` – for example the measurement table
  of a `GSM` or the probe annotation of a `GPL`.

- `Columns(object)`:

  A `data.frame` describing the columns of `Table(object)`.

- `dataTable(object)`:

  The underlying `GEODataTable` object, which holds both `Table()` and
  `Columns()`.

- `GSMList(object)`:

  For a `GSE`, the list of its `GSM` sample objects.

- `GPLList(object)`:

  For a `GSE`, the list of its `GPL` platform objects.

## See also

[`GEOData-class`](http://seandavi.github.io/GEOquery/reference/GEOData-class.md),
[`getGEO`](http://seandavi.github.io/GEOquery/reference/getGEO.md)

## Author

Sean Davis

## Examples

``` r
if (FALSE) { # \dontrun{
  gsm <- getGEO("GSM11805")
  Meta(gsm)$title
  head(Table(gsm))
  Columns(gsm)

  gse <- getGEO("GSE781", GSEMatrix = FALSE)
  names(GSMList(gse))
} # }
```

# Coerce a GEOquery ExpressionSet to a SummarizedExperiment

A thin wrapper around
[`SummarizedExperiment::makeSummarizedExperimentFromExpressionSet()`](https://rdrr.io/pkg/SummarizedExperiment/man/makeSummarizedExperimentFromExpressionSet.html)
used by `getGEO(..., returnType = "SummarizedExperiment")`, and
available directly so existing ExpressionSet results can be modernized
without re-downloading.

## Usage

``` r
as_SummarizedExperiment(eset)
```

## Arguments

- eset:

  An `ExpressionSet`, e.g. an element returned by
  [`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
  for a GSE Series Matrix file.

## Value

A `SummarizedExperiment`.

## Examples

``` r
if (FALSE) { # \dontrun{
  gse <- getGEO("GSE2553")[[1]]
  se <- as_SummarizedExperiment(gse)
} # }
```

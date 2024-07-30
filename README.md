## Status

<!-- badges: start -->
[![R-CMD-check](https://github.com/seandavi/GEOquery/workflows/R-CMD-check/badge.svg)](https://github.com/seandavi/GEOquery/actions)
[![R-CMD-check](https://github.com/seandavi/GEOquery/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seandavi/GEOquery/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Installation

To install from Bioconductor, use the following code:

```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("GEOquery")
```

To install directly from github:

```{r}
library(devtools)
install_github('seandavi/GEOquery')
```

## Usage

See the full vignette in [rmarkdown](https://github.com/seandavi/GEOquery/blob/master/vignettes/GEOquery.Rmd) or visit Bioconductor for details:

- [Release version](http://www.bioconductor.org/packages/release/bioc/html/GEOquery.html)
- [Devel version](http://www.bioconductor.org/packages/devel/bioc/html/GEOquery.html)

## How to contribute

Contributions to GEOquery development can be submitted as a [pull request](https://github.com/seandavi/GEOquery/pulls) or a [feature request issue](https://github.com/seandavi/GEOquery/issues). We recommend following the [Bioconductor coding standards](https://contributions.bioconductor.org/r-code.html) where possible.  

|     Servicd     |   Status      |
|-----------------|---------------|
| Travis (Linux)         | [![Travis Build Status](https://api.travis-ci.org/seandavi/GEOquery.svg?branch=master)](https://travis-ci.org/seandavi/GEOquery) |
| CodeCov         | [![codecov](https://codecov.io/gh/seandavi/GEOquery/branch/master/graph/badge.svg)](https://codecov.io/gh/seandavi/GEOquery) |
| Appveyor (Windows)        | [![AppVeyor](https://ci.appveyor.com/api/projects/status/32r7s2skrgm9ubva?svg=true)](https://ci.appveyor.com/api/projects/status/32r7s2skrgm9ubva?svg=true) |




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
install_github('GEOquery','seandavi')
```

## Usage

See the full vignette in [rmarkdown](https://github.com/seandavi/GEOquery/blob/master/vignettes/GEOquery.Rmd) or visit Bioconductor for details:

- [Release version](http://www.bioconductor.org/packages/release/bioc/html/GEOquery.html)
- [Devel version](http://www.bioconductor.org/packages/devel/bioc/html/GEOquery.html)

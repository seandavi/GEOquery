\name{GEODataTable-class}
\docType{class}
\alias{GEODataTable-class}
\alias{Columns,GEODataTable-method}
\alias{Table,GEODataTable-method}
\alias{show,GEODataTable-method}

\title{Class "GEODataTable" }
\description{Class for storage of column information and data table from
GEO objects}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("GEODataTable",
    ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{columns}:}{Object of class \code{"data.frame"} that
      holds the column descriptions, one row for each column in the table}
    \item{\code{table}:}{Object of class \code{"data.frame"} that
      contains the actual data from a GEO data table (i.e., numbers)}
  }
}
\section{Methods}{
  \describe{
    \item{Columns}{\code{signature(object = "GEODataTable")}: return the
      columns slot from the GEODataTable (i.e., the column descriptions)}
    \item{Table}{\code{signature(object = "GEODataTable")}: return the
      data part of the data table}
    \item{show}{\code{signature(object = "GEODataTable")}: show method
      for GEODataTables }
  }
}
\author{Sean Davis}

\seealso{
  \code{\link{GEOData-class}}
}
\examples{
}
\keyword{classes}

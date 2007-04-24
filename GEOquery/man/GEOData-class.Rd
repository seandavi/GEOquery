\name{GEOData-class}
\docType{class}
\alias{GEOData-class}
\alias{Accession,GEOData-method}
\alias{Columns,GEOData-method}
\alias{Meta,GEOData-method}
\alias{Table,GEOData-method}
\alias{dataTable,GEOData-method}
\alias{show,GEOData-method}

\title{Class "GEOData"}
\description{A virtual class for holding GEO samples, platforms, and datasets}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GEOData", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{header}:}{Object of class \code{"list"} containing metadata}
  }
}
\section{Methods}{
  \describe{
    \item{Accession}{\code{signature(object = "GEOData")}: returns the
      GEO acccession for the current object}
    \item{Columns}{\code{signature(object = "GEOData")}: returns the
      column descriptions for the current object}
    \item{Meta}{\code{signature(object = "GEOData")}: returns the
      metadata for the current object}
    \item{Table}{\code{signature(object = "GEOData")}: returns the
      "Table" for the current object}
    \item{dataTable}{\code{signature(object = "GEOData")}: returns the
      dataTable (column info and data) for the current object}
    \item{show}{\code{signature(object = "GEOData")}: a convenience
      method for showing a GEO object}
	 }
}
\author{ Sean Davis }
\seealso{
	\code{\link{GDS-class}},	\code{\link{GPL-class}},	\code{\link{GSM-class}},	\code{\link{GEODataTable-class}},
}
\keyword{classes}

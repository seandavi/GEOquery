\name{GEODataTable-class}
\docType{class}
\alias{GEODataTable-class}
\alias{Accession,GEODataTable-method}
\alias{Columns,GEODataTable-method}
\alias{Meta,GEODataTable-method}
\alias{Table,GEODataTable-method}
\alias{dataTable,GEODataTable-method}
\alias{show,GEODataTable-method}

\title{Class "GEODataTable"}
\description{Contains the column descriptions and data for the datatable
  part of a GEO object}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GEODataTable", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{columns}:}{Object of class \code{"data.frame"}
      containing information about the columns in the datatable}
    \item{\code{table}:}{Object of class \code{"data.frame"} containing
      the actual data}
  }
}
\section{Methods}{
  \describe{
    \item{Columns}{\code{signature(object = "GEODataTable")}: get the
      column portion of the GEODataTable}
    \item{Table}{\code{signature(object = "GEODataTable")}: get the
      table portion of the GEODataTable}
    \item{show}{\code{signature(object = "GEODataTable")}: convenience
      show method}
	 }
}
\author{ Sean Davis}
\keyword{classes}

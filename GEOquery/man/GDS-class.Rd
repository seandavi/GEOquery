\name{GDS-class}
\docType{class}
\alias{GDS-class}

\title{Class "GDS" }
\description{A class describing a GEO GDS entity}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("GDS", ...)}
  }
\section{Slots}{
	 \describe{
    \item{\code{gpl}:}{Object of class \code{"GPL"}}
    \item{\code{dataTable}:}{Object of class \code{"GEODataTable"}
      containing the data table for the GDS}
    \item{\code{header}:}{Object of class \code{"list"} containing the
      metadata for the GDS; can be accessed via the \code{Meta} accessor}
  }
}
\section{Extends}{
Class \code{"GEOData"}, directly.
}
\section{Methods}{
No methods defined with class "GDS" in the signature, but methods
applying to GEOData are also applicable here.
}
\author{ Sean Davis }

\seealso{
	\code{\link{GEOData-class}}
}
\keyword{classes}

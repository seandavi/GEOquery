\name{GEOData-class}
\docType{class}
\alias{GEOData-class}
\alias{Accession,GEOData-method}
\alias{Columns,GEOData-method}
\alias{Meta,GEOData-method}
\alias{Table,GEOData-method}
\alias{Table<-,GEOData-method}
\alias{dataTable,GEOData-method}
\alias{show,GEOData-method}

\title{Class "GEOData" }
\description{ A virtual class for describing GEO data }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GEOData", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{header}:}{Object of class \code{"list"} }
  }
}
\section{Methods}{
  \describe{
    \item{Accession}{\code{signature(object = "GEOData")}: ... }
    \item{Columns}{\code{signature(object = "GEOData")}: ... }
    \item{Meta}{\code{signature(object = "GEOData")}: ... }
    \item{Table}{\code{signature(object = "GEOData")}: ... }
    \item{dataTable}{\code{signature(object = "GEOData")}: ... }
    \item{show}{\code{signature(object = "GEOData")}: ... }
  }
}
\author{Sean Davis}
\seealso{
  \code{\link{getGEO}}
}
\examples{
gds <- getGEO('GDS505')
Accession(gds)
Meta(gds)
Columns(gds)
gds
}
\keyword{classes}

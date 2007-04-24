\name{GSE-class}
\docType{class}
\alias{GSE-class}
\alias{GPLList,GSE-method}
\alias{GSMList,GSE-method}
\alias{Meta,GSE-method}

\title{Class "GSE"}
\description{Contains a GEO Series entity}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GSE", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{header}:}{Object of class \code{"list"} containing
      metadata for the series}
    \item{\code{gsms}:}{Object of class \code{"list"} containing a list
      of items of class \code{GSM} associated with the series}
    \item{\code{gpls}:}{Object of class \code{"list"} containing a list
      of items of class \code{GPL} associate with the series}
  }
}
\section{Methods}{
  \describe{
    \item{GPLList}{\code{signature(object = "GSE")}: returns a list,
      each item of the list being a \code{GPL} object}
    \item{GSMList}{\code{signature(object = "GSE")}: returns a list,
      each item of the list being a \code{GSM} object}
    \item{Meta}{\code{signature(object = "GSE")}: returns a list, the
      metadata associated with the GSE }
	 }
}
\author{ Sean Davis}
\seealso{
	\code{\link{GPL-class}},\code{\link{GSM-class}}
}
\keyword{classes}

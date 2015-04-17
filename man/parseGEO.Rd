\name{parseGEO}
\alias{parseGEO}
\alias{parseGPL}
\alias{parseGSE}
\alias{parseGDS}
\alias{parseGSM}
\title{Parse GEO text}
\description{
  Workhorse GEO parsers.
}
\usage{
parseGEO(fname, GSElimits, destdir=tempdir(), AnnotGPL=FALSE, getGPL=TRUE)
parseGPL(fname)
parseGDS(fname)
parseGSE(fname, GSElimits)
parseGSM(fname)
}
\arguments{
  \item{fname}{The filename of a SOFT format file.  If the filename ends in .gz, 
  a gzfile() connection is used to read the file directly.}
  \item{GSElimits}{Used to limit the number of GSMs parsed into the GSE
    object; useful for memory management for large GSEs.}
  \item{destdir}{The destination directory into which files will be saved (to be used for caching)}
  \item{AnnotGPL}{Fetch the annotation GPL if available}
  \item{getGPL}{Fetch the GPL associated with a GSEMatrix entity (should remain TRUE for all normal use cases)}
}
\details{
  These are probably not useful to the end-user.  Use getGEO to access
  these functions.  parseGEO simply delegates to the appropriate specific parser.  There
  should be no reason to use the parseGPL, parseGDS, parseGSE, or
  parseGSM functions directly.
}
\value{
  parseGEO returns an object of the associated type.  For example, if it
  is passed the text from a GDS entry, a GDS object is returned.
}
\author{Sean Davis}
\seealso{\code{\link{getGEO}}}

\keyword{ IO }

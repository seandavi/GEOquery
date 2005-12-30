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
parseGEO(con)
parseGPL(txt)
parseGDS(txt)
parseGSE(con)
parseGSM(txt)
}
\arguments{
  \item{txt}{The text from a GEO SOFT format file}
  \item{con}{A connection that contains the body of a GEO SOFT format file}
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

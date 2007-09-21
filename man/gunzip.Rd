\name{gunzip}
\alias{gunzip}
\title{Gunzip a file}
\description{
  gunzip a file
}
\usage{
gunzip(filename, destname = gsub("[.]gz$", "", filename), overwrite = FALSE, remove = TRUE, BFR.SIZE = 1e+07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{The filename to be unzipped}
  \item{destname}{The destination file}
  \item{overwrite}{Boolean indicating whether or not to overwrite a
    destfile of the same name}
  \item{remove}{Boolean indicating whether or not to remove the original
  file after completion}
  \item{BFR.SIZE}{The size of the read buffer....}
}
\details{
  This function was stripped out of R.utils due to breaking some stuff
  on the bioconductor build machine.
}
\value{
  Invisibly, the number of bytes read.
}
\author{Original author: Henrik Bengtsson}
\seealso{\code{\link{gzfile}}}
\keyword{IO}

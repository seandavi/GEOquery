\name{getGEOSuppFiles}
\alias{getGEOSuppFiles}
\title{Get Supplemental Files from GEO}
\description{
  NCBI GEO allows supplemental files to be attached to GEO Series (GSE), GEO
  platforms (GPL), and GEO samples (GSM).  This function "knows" how to
  get these files based on the GEO accession.  No parsing of the
  downloaded files is attempted, since the file format is not generally
  knowable by the computer.
}
\usage{
getGEOSuppFiles(GEO, makeDirectory = TRUE, baseDir = getwd())
}
\arguments{
  \item{GEO}{A GEO accession number such as GPL1073 or GSM1137}
  \item{makeDirectory}{Should a "subdirectory" for the downloaded files
    be created?  Default is TRUE.  If FALSE, the files will be
    downloaded directly into the baseDir.}
  \item{baseDir}{The base directory for the downloads.  Default is the
    current working directory.}
}
\details{
  Again, just a note that the files are simply downloaded.  
}
\value{
  A data frame is returned invisibly with rownames representing the full
  path of the resulting downloaded files and the records in the
  data.frame the output of file.info for each downloaded file.
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\examples{
# a <- getGEOSuppFiles('GSM1137')
# a
}
\keyword{ IO }
\keyword{ database }
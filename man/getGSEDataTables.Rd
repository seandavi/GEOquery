\name{getGSEDataTables}
\alias{getGSEDataTables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get GSE data tables from GEO into R data structures.
}
\description{
In some cases, instead of individual sample records (GSM) containing information regarding sample phenotypes, the GEO Series contains that information in an attached data table.  And example is given by GSE3494 where there are two data tables with important information contained within them.  Using getGEO with the standard parameters downloads the GSEMatrix file which, unfortunately, does not contain the information in the data tables.  This function simply downloads the ``header'' information from the GSE record and parses out the data tables into R data.frames.  
}
\usage{
getGSEDataTables(GSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GSE}{
The GSE identifier, such as ``GSE3494''.
}
}
\value{
A list of data.frames.
}
\references{
http://www.ncbi.nlm.nih.gov/geo
}
\author{
Sean Davis <sdavis2@mail.nih.gov>
}
\seealso{
\code{\link{getGEO}}
}
\examples{
dfl = getGSEDataTables("GSE3494")
lapply(dfl,head)
}
\keyword{ IO }

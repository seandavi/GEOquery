\name{queryGEO}
\alias{queryGEO}
\title{
Query NCBI GEO metadata
}
\description{
This function queries NCBI GEO metadata using the EUtils functionality offered by the \code{\link[rentrez]{entrez_search}} function and returns a \code{data.frame} of results.  
}
\usage{
queryGEO(term, retmax=500, ...)
}
\arguments{
  \item{term}{
  A single string representing a search as described \href{http://www.ncbi.nlm.nih.gov/geo/info/qqtutorial.html}{here}.
}
  \item{retmax}{
  limits the number of returned records.
}
  \item{\dots}{
  passed directly to \code{\link[rentrez]{entrez_search}}
}
}
\details{
This function searches the NCBI GEO \code{gds} database. The details of the \code{term} argument are described \href{http://www.ncbi.nlm.nih.gov/geo/info/qqtutorial.html}{on the GEO documentation site}. The function will return all query results which may include GDS, GPL, GSE, or GSM records. The retmax parameter will limit the number of results, so you may need to adjust this value higher if your query returns a lot of results.
}
\value{
A data frame with the following columns:
  \item{accession}{The NCBI GEO accession, appropriate for use in a call to \code{\link{getGEO}}}
}
\references{
\url{http://www.ncbi.nlm.nih.gov/geo/}
}
\author{
Sean Davis
}
\seealso{
\code{\link{getGEO}}, \code{\link[rentrez]{entrez_search}}
}
\examples{
res = queryGEO('sarcoma AND "cel"[Supplementary Files] AND 20:20000[Number of Samples]')
dim(res)
head(res)
}


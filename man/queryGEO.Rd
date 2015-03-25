\name{queryGEO}
\alias{queryGEO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Query NCBI GEO metadata
}
\description{
This function queries NCBI GEO \textit{metadata} using the EUtils functionality offered by the \code{\link[rentrez]{entrez_search}} function and returns a \code{data.frame} of results.  
}
\usage{
queryGEO(term, retmax=500, ...)
}
%- maybe also 'usage' for other objects documented here.
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
This function searches the NCBI GEO \code{gds} database. The details of the \code{term} argument are described \href{http://www.ncbi.nlm.nih.gov/geo/info/qqtutorial.html}{on the GEO documentation site}. The function will return all query results which may include GDS, GPL, GSE, or GSM records. The retmax parameter will \textit{limit} the number of results, so you may need to adjust this value higher if your query returns a lot of results.
}
\value{
A data frame with the following columns:
  \item{accession}{The NCBI GEO accession, appropriate for use in a call to \code{\link{getGEO}}}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (term, retmax=500, ...) 
{
    require(rentrez)
    tmp = entrez_summary(db = "gds", id = entrez_search(db = "gds", 
        term, ...)$ids)
    return(data.frame(do.call(rbind, lapply(tmp, function(y) {
        return(c(accession = y$Accession, title = y$title, taxon = y$taxon, 
            entityType = y$entryType, pubDate = y$PDAT, nSamples = y$n_samples))
    }))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

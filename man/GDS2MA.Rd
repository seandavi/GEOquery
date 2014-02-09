\name{Converting}
\alias{GDS2MA}
\alias{GDS2eSet}
\title{Convert a GDS data structure to a BioConductor data structure}
\description{
  Functions to take a GDS data structure from getGEO and coerce it to
  limma MALists or ExpressionSets.  
}
\usage{
GDS2MA(GDS,do.log2=FALSE,GPL=NULL,AnnotGPL=TRUE,getGPL=TRUE)
GDS2eSet(GDS,do.log2=FALSE,GPL=NULL,AnnotGPL=TRUE,getGPL=TRUE)
}
\arguments{
  \item{GDS}{The GDS datastructure returned by getGEO}
  \item{do.log2}{Boolean, should the data in the GDS be log2 transformed
    before inserting into the new data structure}
  \item{GPL}{Either a GPL data structure (from a call to getGEO) or
    NULL.  If NULL, this will cause a call to getGEO to produce a GPL.
    The gene information from the GPL is then used to construct the
    \code{genes} slot of the resulting limma \code{MAList} object or the
    \code{featureData} slot of the \code{ExpressionSet} instance. }
  \item{AnnotGPL}{In general, the annotation GPL files will be available
    for GDS records, so the default is to use these files over the
    user-submitted GPL files}
  \item{getGPL}{A boolean defaulting to TRUE as to whether or not to
    download and include GPL information when converting to
    ExpressionSet or MAList.
    You may want to set this to FALSE if you know that you are going to
    annotate your featureData using Bioconductor tools rather than relying
    on information provided through NCBI GEO.  Download times can also
    be greatly reduced by specifying FALSE.}
}
\details{
  This function just rearranges one data structure into another.  For
  GDS, it also deals appropriately with making the "targets" list item
  for the limma data structure and the phenoData slot of ExpressionSets.
}
\value{
  \item{GDS2MA}{A limma MAList}
  \item{GDS2eSet}{An ExpressionSet object}
}
\references{See the limma and ExpressionSet help in the appropriate packages}
\author{Sean Davis}
\examples{

\dontrun{gds505 <- getGEO('GDS505')}
\dontrun{MA <- GDS2MA(gds505)}
\dontrun{eset <- GDS2eSet(gds505)}

}
\keyword{IO}

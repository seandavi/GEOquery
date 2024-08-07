% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq.R
\name{readRNAQuantAnnotation}
\alias{readRNAQuantAnnotation}
\title{Read RNA-seq quantification annotation from GEO}
\usage{
readRNAQuantAnnotation(link)
}
\arguments{
\item{link}{A link to the annotation file}
}
\value{
A data frame of annotation information with gene IDs as row names
}
\description{
This function reads the annotation file from a GEO link. The annotation
file is expected to be a tab-separated file with the first column
containing the gene IDs and the remaining columns containing the
annotation information.
}
\keyword{internal}

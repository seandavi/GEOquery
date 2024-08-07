% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq.R
\name{getRNASeqQuantGenomeInfo}
\alias{getRNASeqQuantGenomeInfo}
\title{Extract genome build and species for GEO RNA-seq quantification}
\usage{
getRNASeqQuantGenomeInfo(gse)
}
\arguments{
\item{gse}{GEO accession number}
}
\value{
A character vector with the genome build and species information
}
\description{
This function extracts the genome build and species information
for a GEO RNA-seq quantification.
}
\examples{
extractGenomeBuildSpecies("GSE164073")

}

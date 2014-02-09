\name{getGEO}
\alias{getGEO}
\title{Get a GEO object from NCBI or file}
\description{
  This function is the main user-level function in the GEOquery
  package.  It directs the download (if no filename is specified) and
  parsing of a GEO SOFT format file into an R data structure
  specifically designed to make access to each of the important parts of
  the GEO SOFT format easily accessible.
}
\usage{
getGEO(GEO = NULL, filename = NULL, destdir = tempdir(), GSElimits=NULL,
GSEMatrix=TRUE,AnnotGPL=FALSE,getGPL=TRUE)
}
\arguments{
  \item{GEO}{A character string representing a GEO object for download
    and parsing.  (eg., 'GDS505','GSE2','GSM2','GPL96')}
  \item{filename}{The filename of a previously downloaded GEO SOFT format file or its gzipped
    representation (in which case the filename must end in .gz).  Either
  one of GEO or filename may be specified, not both.  GEO series matrix
  files are also handled.  Note that since a single file is being parsed,
  the return value is not a list of esets, but a single eset when GSE matrix
  files are parsed.}
  \item{destdir}{The destination directory for any downloads.  Defaults
    to the architecture-dependent tempdir.  You may want to specify a
    different directory if you want to save the file for later use.
    Doing so is a good idea if you have a slow connection, as some of
    the GEO files are HUGE!}
  \item{GSElimits}{This argument can be used to load only a contiguous
    subset of the GSMs from a GSE.  It should be specified as a vector
    of length 2 specifying the start and end (inclusive) GSMs to load.
    This could be useful for splitting up large GSEs into more
    manageable parts, for example.}
  \item{GSEMatrix}{A boolean telling GEOquery whether or not to use GSE
    Series Matrix files from GEO.  The parsing of these files can be
    many orders-of-magnitude faster than parsing the GSE SOFT format
    files.  Defaults to TRUE, meaning that the SOFT format parsing will
    not occur; set to FALSE if you for some reason need other columns
    from the GSE records.}
  \item{AnnotGPL}{A boolean defaulting to FALSE as to whether or not to
    use the Annotation GPL information.  These files are nice to use
    because they contain up-to-date information remapped from Entrez
    Gene on a regular basis.  However, they do not exist for all GPLs;
    in general, they are only available for GPLs referenced by a GDS}
  \item{getGPL}{A boolean defaulting to TRUE as to whether or not to
    download and include GPL information when getting a GSEMatrix file.
    You may want to set this to FALSE if you know that you are going to
    annotate your featureData using Bioconductor tools rather than relying
    on information provided through NCBI GEO.  Download times can also
    be greatly reduced by specifying FALSE.}
}
\details{
  getGEO functions to download and parse information available from NCBI
  GEO (\url{http://www.ncbi.nlm.nih.gov/geo}).  Here are some details
  about what is avaible from GEO.  All entity types are handled by
  getGEO and essentially any information in the GEO SOFT format is
  reflected in the resulting data structure.
  
  From the GEO website:

  The Gene Expression Omnibus (GEO) from NCBI serves as a public
  repository for a wide range of high-throughput experimental
  data. These data include single and dual channel microarray-based
  experiments measuring mRNA, genomic DNA, and protein abundance, as
  well as non-array techniques such as serial analysis of gene
  expression (SAGE), and mass spectrometry proteomic data. At the most
  basic level of organization of GEO, there are three entity types that
  may be supplied by users: Platforms, Samples, and Series.
  Additionally, there is a curated entity called a GEO dataset.

  A Platform record describes the list of elements on the array (e.g.,
  cDNAs, oligonucleotide probesets, ORFs, antibodies) or the list of
  elements that may be detected and quantified in that experiment (e.g.,
  SAGE tags, peptides). Each Platform record is assigned a unique and
  stable GEO accession number (GPLxxx). A Platform may reference many
  Samples that have been submitted by multiple submitters. 

  A Sample record describes the conditions under which an individual
  Sample was handled, the manipulations it underwent, and the abundance
  measurement of each element derived from it. Each Sample record is
  assigned a unique and stable GEO accession number (GSMxxx). A Sample
  entity must reference only one Platform and may be included in
  multiple Series.
  
  A Series record defines a set of related Samples considered to be part
  of a group, how the Samples are related, and if and how they are
  ordered. A Series provides a focal point and description of the
  experiment as a whole. Series records may also contain tables
  describing extracted data, summary conclusions, or analyses. Each
  Series record is assigned a unique and stable GEO accession number
  (GSExxx). 

  GEO DataSets (GDSxxx) are curated sets of GEO Sample data. A GDS
  record represents a collection of biologically and statistically
  comparable GEO Samples and forms the basis of GEO's suite of data
  display and analysis tools. Samples within a GDS refer to the same
  Platform, that is, they share a common set of probe elements. Value
  measurements for each Sample within a GDS are assumed to be calculated
  in an equivalent manner, that is, considerations such as background
  processing and normalization are consistent across the
  dataset. Information reflecting experimental design is provided
  through GDS subsets.

}
\value{
  An object of the appropriate class (GDS, GPL, GSM, or GSE) is
  returned.  If the GSEMatrix option is used, then a list of
  ExpressionSet objects is returned, one for each SeriesMatrix file
  associated with the GSE accesion.  If the filename argument is used
  in combination with a GSEMatrix file, then the return value is a single
  ExpressionSet.
}
\author{Sean Davis}
\section{Warning }{Some of the files that are downloaded, particularly
  those associated with GSE entries from GEO are absolutely ENORMOUS and
  parsing them can take quite some time and memory.  So, particularly
  when working with large GSE entries, expect that you may need a good
  chunk of memory and that coffee may be involved when parsing....}

\seealso{\code{\link{getGEOfile}}}
\examples{
# gds <- getGEO("GDS10")
# gds
}

\keyword{IO}

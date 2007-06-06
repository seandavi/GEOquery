getGEO <- function(GEO=NULL,
                   filename=NULL,
                   destdir=tempdir(),
                   GSElimits=NULL,GSEMatrix=FALSE) {
  filename <- filename
  if(!is.null(GSElimits)) {
    if(length(GSElimits)!=2) {
      stop('GSElimits should be an integer vector of length 2, like (1,10) to include GSMs 1 through 10')
    }
  }
  if(is.null(GEO) & is.null(filename)) {
    stop("You must supply either a filename of a GEO file or a GEO accession")
  }
  if(is.null(filename)) {
    GEO <- toupper(GEO)
    if(GSEMatrix) {
      return(getAndParseGSEMatrices(GEO))
    }
    filename <- getGEOfile(GEO,destdir=destdir)
  }
  if(!is.null(filename) & GSEMatrix) {
    stop("Currently, getting GSEmatrix from local files is not supported")
  }
  if(length(grep('\\.gz$',filename,perl=TRUE))>0) {
    con <- gzfile(filename,'r')
  } else {
    con <- file(filename,'r')
  }
  ret <- parseGEO(con,GSElimits)
  close(con)
  return(ret)
}       

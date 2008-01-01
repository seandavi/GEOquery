getGEO <- function(GEO=NULL,
                   filename=NULL,
                   destdir=tempdir(),
                   GSElimits=NULL,GSEMatrix=TRUE) {
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
### I was having problems with gunzip on Windows, which is why
###  I added the lines below.  I'm giving it another try to see if
###  something has changed with gunzip.  
#      if(.Platform$OS.type=='windows') {
#        warning('GSEMatrix parsing not currently supported on Windows')
#        return(NULL)
#      }
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

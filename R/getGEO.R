getGEO <- function(GEO=NULL,
                   filename=NULL,
                   destdir=tempdir(),
                   GSElimits=NULL,GSEMatrix=TRUE,
                   AnnotGPL=FALSE,
                   getGPL=TRUE) {
  con <- NULL
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
    geotype <- toupper(substr(GEO,1,3))
    if(GSEMatrix & geotype=='GSE') {
      return(getAndParseGSEMatrices(GEO,destdir,AnnotGPL=AnnotGPL,getGPL=getGPL))
    }
    filename <- getGEOfile(GEO,destdir=destdir,AnnotGPL=AnnotGPL)
  }      
  ret <- parseGEO(filename,GSElimits,destdir,AnnotGPL=AnnotGPL,getGPL=getGPL)
  return(ret)
}       

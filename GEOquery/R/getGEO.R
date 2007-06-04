getGEO <- function(GEO=NULL,
                   filename=NULL,
                   destdir=tempdir(),
                   GSElimits=NULL) {
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
    filename <- getGEOfile(GEO,destdir=destdir)
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
                   
                   

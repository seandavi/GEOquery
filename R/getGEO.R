getGEO <- function(GEO=NULL,
                   filename=NULL,
                   destdir=tempdir()) {
  filename <- filename
  if(is.null(GEO) & is.null(filename)) {
    stop("You must supply either a filename of a GEO file or a GEO accession")
  }
  if(is.null(filename)) {
    GEO <- toupper(GEO)
    filename <- getGEOfile(GEO,destdir=destdir)
  }
  writeLines(filename)
  if(length(grep('\\.gz$',filename,perl=TRUE))>0) {
    con <- gzfile(filename,'r')
  } else {
    con <- file(filename,'r')
  }
  return(parseGEO(con))
  close(con)
}
                   
                   

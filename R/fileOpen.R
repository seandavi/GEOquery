fileOpen <- function(fname,...) {
  con <- NULL
  if(length(grep('\\.gz$',fname,perl=TRUE))>0) {
    con <- gzfile(fname,open='rt')
  } else {
    con <- file(fname,'r')
  }
  return(con)
}

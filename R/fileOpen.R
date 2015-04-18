fileOpen <- function(fname,...) {
  con <- NULL
  if(!file.exists(fname)) {
    stop(sprintf("File %s does not appear to exist.",fname))
  }
  if(length(grep('\\.gz$',fname,perl=TRUE))>0) {
    con <- gzfile(fname,open='rt')
  } else {
    con <- file(fname,'r')
  }
  return(con)
}

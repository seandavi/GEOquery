getDirListing <- function(url) {
  print(url)
  a <- getURL(url)
  b <- textConnection(a)
  d <- read.table(b,header=FALSE)
  close(b)
  return(d)
}

getGEOSuppFiles <- function(GEO,makeDirectory=TRUE,baseDir=getwd()) {
  geotype <- toupper(substr(GEO,1,3))
  storedir <- baseDir
  if(makeDirectory) {
    try(dir.create(GEO))
    storedir <- file.path(baseDir,GEO)
  }
  fileinfo <- list()
  if(geotype=='GSM') {
    gsmstub <- gsub('\\d{1,3}$','nnn',GEO,perl=TRUE)
    url <- sprintf("ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/%s/%s/",gsmstub,GEO)
  }
  if(geotype=='GSE') {
    url <- sprintf("ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/series/%s/",GEO)
  }
  if(geotype=='GPL') {
    url <- sprintf("ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/platforms/%s/",GEO)
  }
  dirlist <- getDirListing(url)
  fnames <- as.character(dirlist[,9])
  for(i in fnames) {
    download.file(file.path(url,i),
                  destfile=file.path(storedir,i),
                  method=getOption('download.file.method.GEOquery'))
    fileinfo[[file.path(storedir,i)]] <- file.info(file.path(storedir,i))
  }
  invisible(do.call(rbind,fileinfo))
}
    

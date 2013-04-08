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
  fileinfo <- list()
  stub = gsub('\\d{1,3}$','nnn',GEO,perl=TRUE)
  if(geotype=='GSM') {
    url <- sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/",stub,GEO)
  }
  if(geotype=='GSE') {
    url <- sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/",stub,GEO)
  }
  if(geotype=='GPL') {
    url <- sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/platform/%s/%s/suppl/",stub,GEO)
  }
  dirlist <- try(getDirListing(url),silent=TRUE)
  if(inherits(dirlist,'try-error')) {
    message('No supplemental files found')
    return(NULL)
  }
  if(makeDirectory) {
    try(dir.create(GEO),silent=TRUE)
    storedir <- file.path(baseDir,GEO)
  }
  fnames <- as.character(dirlist[,9])
  for(i in fnames) {
    download.file(file.path(url,i),
                  destfile=file.path(storedir,i),
                  mode='wb',
                  method=getOption('download.file.method.GEOquery'))
    fileinfo[[file.path(storedir,i)]] <- file.info(file.path(storedir,i))
  }
  invisible(do.call(rbind,fileinfo))
}
    

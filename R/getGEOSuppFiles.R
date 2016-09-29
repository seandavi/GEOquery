getDirListing <- function(url) {
  message(url)
  # Takes a URL and returns a character vector of filenames
  a <- getURL(url)
  ## Renaud Gaujoux reported problems behind firewall
  ## where the ftp index was converted to html content
  ## The IF statement here is his fix--harmless for the rest
  ## of us.
  if( grepl("<HTML", a, ignore.case=T) ){ # process HTML content
    doc <- XML::htmlParse(a)
    links <- XML::xpathSApply(doc, "//a/@href")
    XML::free(doc)
    
    b <- as.matrix(links)
    message('OK')
  } else { # standard processing of txt content
    tmpcon <- textConnection(a, "r")
    b <- read.table(tmpcon)
    close(tmpcon)
  }
  b <- as.character(b[,ncol(b)])
  return(b)
}

getGEOSuppFiles <- function(GEO,makeDirectory=TRUE,baseDir=getwd()) {
  geotype <- toupper(substr(GEO,1,3))
  storedir <- baseDir
  fileinfo <- list()
  stub = gsub('\\d{1,3}$','nnn',GEO,perl=TRUE)
  if(geotype=='GSM') {
    url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/",stub,GEO)
  }
  if(geotype=='GSE') {
    url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/",stub,GEO)
  }
  if(geotype=='GPL') {
    url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platform/%s/%s/suppl/",stub,GEO)
  }
  fnames <- try(getDirListing(url),silent=TRUE)
  if(inherits(fnames,'try-error')) {
    message('No supplemental files found.')
    message('Check URL manually if in doubt')
    message(url)
    return(NULL)
  }
  if(makeDirectory) {
    suppressWarnings(dir.create(storedir <- file.path(baseDir,GEO)))
  }
  for(i in fnames) {
    download.file(file.path(url,i),
                  destfile=file.path(storedir,i),
                  mode='wb',
                  method=getOption('download.file.method.GEOquery'))
    fileinfo[[file.path(storedir,i)]] <- file.info(file.path(storedir,i))
  }
  invisible(do.call(rbind,fileinfo))
}
    

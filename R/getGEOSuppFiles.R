#' get a directory listing from NCBI GEO
#'
#' This one makes some assumptions about the
#' structure of the HTML response returned.
#'
#' @importFrom xml2 xml_find_all, xml_text, read_html
#' 
getDirListing <- function(url) {
  # Takes a URL and returns a character vector of filenames
    a <- xml2::read_html(url)
    fnames = grep('^G',xml_text(xml_find_all(a,'//a/@href')),value=TRUE)
  return(fnames)
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
    download.file(paste(file.path(url,i),'tool=geoquery',sep="?"),
                  destfile=file.path(storedir,i),
                  mode='wb',
                  method=getOption('download.file.method.GEOquery'))
    fileinfo[[file.path(storedir,i)]] <- file.info(file.path(storedir,i))
  }
  invisible(do.call(rbind,fileinfo))
}
    

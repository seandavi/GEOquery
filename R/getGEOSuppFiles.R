#' get a directory listing from NCBI GEO
#' 
#' This one makes some assumptions about the structure of the HTML response
#' returned.
#' 
#' @importFrom xml2 read_html xml_text xml_find_all 
getDirListing <- function(url) {
  # Takes a URL and returns a character vector of filenames
    a <- xml2::read_html(url)
    fnames = grep('^G',xml_text(xml_find_all(a,'//a/@href')),value=TRUE)
  return(fnames)
}



#' Get Supplemental Files from GEO
#' 
#' NCBI GEO allows supplemental files to be attached to GEO Series (GSE), GEO
#' platforms (GPL), and GEO samples (GSM).  This function "knows" how to get
#' these files based on the GEO accession.  No parsing of the downloaded files
#' is attempted, since the file format is not generally knowable by the
#' computer.
#' 
#' Again, just a note that the files are simply downloaded.
#' 
#' @param GEO A GEO accession number such as GPL1073 or GSM1137
#' @param makeDirectory Should a "subdirectory" for the downloaded files be
#' created?  Default is TRUE.  If FALSE, the files will be downloaded directly
#' into the baseDir.
#' @param baseDir The base directory for the downloads.  Default is the current
#' working directory.
#' @return A data frame is returned invisibly with rownames representing the
#' full path of the resulting downloaded files and the records in the
#' data.frame the output of file.info for each downloaded file.
#' @author Sean Davis <sdavis2@@mail.nih.gov>
#' @keywords IO database
#' @examples
#' 
#' # a <- getGEOSuppFiles('GSM1137')
#' # a
#' 
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
    

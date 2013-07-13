.onLoad <- function(lib,pkg) {
  if(length(grep('linux',R.Version()$os))==1) {
    packageStartupMessage("Setting options('download.file.method.GEOquery'='curl')")
    options('download.file.method.GEOquery'='curl')
  } else {
    packageStartupMessage("Setting options('download.file.method.GEOquery'='auto')")
    options('download.file.method.GEOquery'='auto')
  }
}

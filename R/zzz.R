.onLoad <- function(lib,pkg) {
  require(methods)
  if(length(grep('linux',R.Version()$os))==1) {
    message("Setting options('download.file.method.GEOquery'='curl')")
    options('download.file.method.GEOquery'='curl')
  } else {
    message("Setting options('download.file.method.GEOquery'='curl')")
    options('download.file.method.GEOquery'='auto')
  }
}

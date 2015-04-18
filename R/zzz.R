.onLoad <- function(lib,pkg) {
  packageStartupMessage("Setting options('download.file.method.GEOquery'='auto')")
  options('download.file.method.GEOquery'='auto')
}

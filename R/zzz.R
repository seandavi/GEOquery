.onLoad <- function(lib,pkg) {
  require(methods)
  message("Setting options('download.file.method'='curl')")
  options('download.file.method'='curl')
}

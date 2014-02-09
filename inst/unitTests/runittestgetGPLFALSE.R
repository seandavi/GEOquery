"testCaseGSEgetGPLFALSE" <-
function() {
  gse = getGEO('GSE2553')[[1]]

  checkTrue(validObject(gse))
  checkEqualsNumeric(0,ncol(gse))
  checkEqualsNumeric(12600,nrow(gse))
}

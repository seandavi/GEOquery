"testCaseGSEgetGPLFALSE" <-
function() {
  gse = getGEO('GSE2553',getGPL=FALSE)[[1]]

  checkTrue(validObject(gse))
  checkEqualsNumeric(0,ncol(fData(gse)))
  checkEqualsNumeric(12600,nrow(fData(gse)))
}


"testCaseGDS2eSetgetGPLFALSE" <-
function() {
  gds = getGEO('GDS10')

  eset = GDS2eSet(gds,getGPL=FALSE)
  checkTrue(validObject(eset))
  checkEqualsNumeric(0,ncol(fData(eset)))
  checkEqualsNumeric(39114,nrow(fData(eset)))
}

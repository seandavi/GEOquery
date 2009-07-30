test.getGEO_GDS <- function() {
  gds <- getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))
  checkEquals(class(gds)[[1]],'GDS',msg="GDS class check")
  checkEquals(dim(Table(gds))[1],22645,msg="GDS row count")
  checkEquals(dim(Table(gds))[2],19,msg="GDS column count")
}

test.getGEO_GSM <- function() {
  gsm <- getGEO(filename=system.file("extdata/GSM11805.txt.gz",package="GEOquery"))
  checkEquals(class(gsm)[[1]],'GSM',msg="GSM class check")
  checkEquals(dim(Table(gsm))[1],22283,msg="GSM row count")
  checkEquals(dim(Table(gsm))[2],3,msg="GSM column count")
  checkEquals(class(Meta(gsm)),"list",msg="GSM Meta() return class test")
  checkEquals(length(Meta(gsm)),28,msg="GSM Meta() method test")
}

test.getGEO_GPL <- function() {
  gpl <- getGEO(filename=system.file("extdata/GPL96.txt.gz",package="GEOquery"))
  checkEquals(class(gpl)[[1]],'GPL',msg="GPL class check")
  checkEquals(dim(Table(gpl))[1],22283,msg="GPL row count")
  checkEquals(dim(Table(gpl))[2],16,msg="GPL column count")
  checkEquals(class(Meta(gpl)),"list",msg="GPL Meta() return class test")
  checkEquals(length(Meta(gpl)),25,msg="GPL Meta() method test")
}




## Test unit 'testSOFTFormatGSE'

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runittestSOFTFormatGSE.R"
		.Log$..File <- ""
		.Log$..Obj <- ""
		.Log$..Tag <- ""
		.Log$..Msg <- ""
		rm(..Test, envir = .Log)
	}
}

.tearDown <-
function () {
	## Specific actions for svUnit: clean up context
	if ("package:svUnit" %in% search()) {
		.Log$..Unit <- ""
		.Log$..File <- ""
		.Log$..Obj <- ""
		.Log$..Tag <- ""
		.Log$..Msg <- ""
		rm(..Test, envir = .Log)
	}
}

"testSOFTFormatGSE" <-
function() {
  gse = getGEO('GSE1563',GSEMatrix=FALSE)

  checkEquals(62,length(GSMList(gse)))
  checkEquals(1,length(GPLList(gse)))
  checkEquals(12625,nrow(Table(GPLList(gse)[[1]])))
  checkEquals(12625,nrow(Table(GSMList(gse)[[1]])))
  sapply(GSMList(gse),function(x) {
      checkTrue(inherits(x,'GSM'),'all GSMList members should be "GSM"s')
      checkEquals(12625,nrow(Table(x)),'all GSMs should have the 12625 rows')})
  checkTrue(inherits(GPLList(gse)[[1]],'GPL'))
}

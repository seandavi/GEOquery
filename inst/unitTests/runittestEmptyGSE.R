## Test unit 'testEmptyGSE'

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runittestEmptyGSE.R"
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

"testEmptyGSE" <-
function() {
    gse = getGEO('GSE11413')

    checkTrue(inherits(gse,'list'),'gse should be a list')
    checkTrue(inherits(gse[[1]],'ExpressionSet'))
    checkEquals(12,nrow(pData(gse[[1]])))
    checkEquals(0,nrow(fData(gse[[1]])))
}

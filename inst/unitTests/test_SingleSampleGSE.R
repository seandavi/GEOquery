## Test unit 'testSingleSampleGSE'

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runittestSingleSampleGSE.R"
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

"testSingleSampleGSE" <-
function() {
    gse = getGEO('GSE11595')

    checkTrue(inherits(gse[[1]],'ExpressionSet'),'result should be an ExpressionSet')
    checkEqualsNumeric(1,ncol(gse[[1]]),'Should be one column for one sample')
}

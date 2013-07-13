## Test unit 'testGDSConversions'

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runittestGDSConversions.R"
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

"testGDSConversions" <-
function() {
    gds = getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))

    eset = GDS2eSet(gds)

    checkTrue(inherits(eset,'ExpressionSet'),'eset should be an ExpressionSet')
    checkEquals(pubMedIds(experimentData(eset)),'14641932','basic experimentData check failed')
    checkEqualsNumeric(nrow(eset),22645,'eset has wrong number of rows')
    checkEqualsNumeric(ncol(eset),17,'eset has wrong number of columns')

    malist = GDS2MA(gds)

    checkTrue(inherits(malist,'MAList'),'malist should be an MAList')
    checkEquals(nrow(malist),22645,'malist has wrong number of rows')
    checkEquals(ncol(malist),17,'malist has wrong number of columns')
}

## Test unit 'testSuppFileSupport'

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runittestSuppFileSupport.R"
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

"testSuppFileSupport" <-
function() {
    fres = getGEOSuppFiles('GSE1000')

    checkEquals(10,ncol(fres))
    checkEquals(2,nrow(fres))

    fres = getGEOSuppFiles('GSM15789')

    checkEquals(10,ncol(fres))
    checkEquals(1,nrow(fres))

    fres = getGEOSuppFiles('GSM15789',baseDir=tempdir())

    checkEquals(10,ncol(fres))
    checkEquals(1,nrow(fres))
}

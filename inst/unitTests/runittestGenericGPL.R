## Test unit 'testGenericGPL'

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runittestGenericGPL.R"
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

"testGenericGPL" <-
function() {
        gpl = getGEO('GPL96')

        checkTrue(inherits(gpl,'GPL'),'gpl is not a GPL object!')
        checkTrue(inherits(gpl,'GEOData'),'gpl does not inherit from GEOData!')
        checkEqualsNumeric(nrow(Table(gpl)),22283,'nrow(gpl) does not match!')
        checkEqualsNumeric(ncol(Table(gpl)),16,'ncol(gpl) does not match!')
        checkTrue(inherits(Meta(gpl),'list'),'Meta(gpl) should be a list')
        checkTrue(inherits(Table(gpl),'data.frame'),'Table(gpl) should be a data.frame')
        checkTrue(inherits(dataTable(gpl),'GEODataTable'),'dataTable(gpl) should be a GEODataTable')
    }

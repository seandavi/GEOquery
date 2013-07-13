## Test unit 'testGenericGSM'

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runittestGenericGSM.R"
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

"testGenericGSM" <-
function() {
    gsm = getGEO('GSM11805')

    checkTrue(inherits(gsm,'GSM'))
    checkTrue(inherits(Meta(gsm),'list'))
    checkTrue(inherits(Table(gsm),'data.frame'))
    checkTrue(inherits(dataTable(gsm),'GEODataTable'))
    checkEquals(Accession(gsm),'GSM11805')
    checkEquals(nrow(Table(gsm)),22283)
    checkEquals(ncol(Table(gsm)),3)
    checkEquals(length(Meta(gsm)),28)
}

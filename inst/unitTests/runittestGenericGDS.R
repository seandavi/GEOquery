## Test unit 'testGenericGDS'

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runittestGenericGDS.R"
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

"testGenericGDS" <-
function() {
    gds = getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))

    checkTrue(inherits(gds,'GDS'),'gds should be a "GDS"')
    checkTrue(inherits(Meta(gds),'list'),'Meta(gds) should be a list')
    checkTrue(inherits(Table(gds),'data.frame'),'Table(gds) should be a data.frame')
    checkTrue(inherits(dataTable(gds),'GEODataTable'),'dataTable(gds) should be a GEODataTable')
    checkEquals(nrow(Table(gds)),22645,'gds has wrong number of rows!')
    checkEquals(ncol(Table(gds)),19,'gds has wrong number of columns!')
    checkEquals(length(Meta(gds)),23,'gds has wrong number of Metadata entries!')
}

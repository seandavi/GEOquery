### $Id$

### Generic GEO classes:

setClass("GEODataTable",
         representation(
                        columns="data.frame",
                        table="data.frame"
                        ),
         prototype=list(
           columns=data.frame(matrix(nr=0,nc=0)),
           table=data.frame(matrix(nr=0,nc=0))
           )
         )


setClass("GEOData",
         representation(
                        header="list"
                        ),
         prototype=list(
           header=list()
         )
         )

setClass("GSE",
         representation(
                        header="list",
                        gsms="list",
                        gpls="list"
                        ),
         prototype=list(
           header=list(),
           gsms=list(),
           gpls=list()
         )
         )

setClass("GPL",
         representation(
                        dataTable='GEODataTable'
                        ),
         prototype=list(
           dataTable=new("GEODataTable")
           ),
         contains="GEOData")
setClass("GSM",
         representation(
                        dataTable='GEODataTable'
                        ),
         prototype=list(
           dataTable=new("GEODataTable")
           ),
         contains="GEOData")
setClass("GDS",
         representation(
                        gpl="GPL",
                        dataTable='GEODataTable'
                        ),
         prototype=list(
           gpl=new("GPL"),
           dataTable=new("GEODataTable")
           ),
         contains="GEOData")

printHead <- function(x)
#  Print leading 5 elements or rows of atomic object
#  From limma and Gordon Smyth
{
        if(is.atomic(x)) {
                d <- dim(x)
                if(length(d)<2) which <- "OneD"
                if(length(d)==2) which <- "TwoD"
                if(length(d)>2) which <- "Array"
        } else {
                if(inherits(x,"data.frame")) {
                        d <- dim(x)
                        which <- "TwoD"
                } else
                        which <- "Recursive"
        }
        switch(which,
        OneD={
                n <- length(x)
                if(n > 20) {
                        print(x[1:5])
                        cat(n-5,"more elements ...\n")
                } else
                        print(x)
        },
        TwoD={
                n <- d[1]
                if(n > 10) {
                        print(x[1:5,])
                        cat(n-5,"more rows ...\n")
                } else
                        print(x)
        },
        Array={
                n <- d[1]
                if(n > 10) {
                        dn <- dimnames(x)
                        dim(x) <- c(d[1],prod(d[-1]))
                        x <- x[1:5,]
                        dim(x) <- c(5,d[-1])
                        if(!is.null(dn[[1]])) dn[[1]] <- dn[[1]][1:5]
                        dimnames(x) <- dn
                        print(x)
                        cat(n-5,"more rows ...\n")
                } else
                        print(x)
        },
        Recursive=print(x)
        )
}

chopColumns <- function(x) {
  apply(x,2,substring,first=1,last=35)
}

setGeneric("Meta",
           function(object) standardGeneric("Meta"))
setGeneric("Accession",
           function(object) standardGeneric("Accession"))
setGeneric("dataTable",
           function(object) standardGeneric("dataTable"))
setGeneric("Columns",
           function(object) standardGeneric("Columns"))
setGeneric("GPL",
           function(object) standardGeneric("GPL"))
setGeneric("GSM",
           function(object) standardGeneric("GSM"))
setGeneric("GPLList",
           function(object) standardGeneric("GPLList"))
setGeneric("GSMList",
           function(object) standardGeneric("GSMList"))
setGeneric("Table",
           function(object) standardGeneric("Table"))

#' @export
setMethod("Meta","GEOData",
          function(object) {
            return(object@header)
          }
          )
#' @export
setMethod("Meta","GSE",
          function(object) {
            return(object@header)
          }
          )
          
#' @export
setMethod("Accession","GEOData",
		 function(object) {
		   return(Meta(object)$geo_accession)
		 }
		 )

#' @export
setMethod("dataTable","GEOData",
          function(object) {
            return(object@dataTable)
          }
          )

#' @export
setMethod("Table","GEODataTable",
          function(object) {
            return(object@table)
          }
          )

#' @export
setMethod("Columns","GEODataTable",
          function(object) {
            return(object@columns)
          }
          )

#' @export
setMethod("Table","GEOData",
          function(object) {
            return(Table(dataTable(object)))
          }
          )

#' @export
setMethod("Columns","GEOData",
          function(object) {
            return(Columns(dataTable(object)))
          }
          )

#' @export
setMethod("GPLList","GSE",
          function(object) {
            return(object@gpls)
          }
          )

#' @export
setMethod("GPL","GDS",
          function(object) {
            return(object@gpl)
          }
          )

#' @export
setMethod("GSMList","GSE",
          function(object) {
            return(object@gsms)
          }
          )
          

#' @export
setMethod("show","GEODataTable",
          function(object) {
            cat("An object of class \"",class(object),"\"\n",sep="")
            cat("****** Column Descriptions ******\n")
            print(Columns(object))
                  cat("****** Data Table ******\n")
            printHead(Table(object))
          })

#' @export
setMethod("show","GEOData",
          function(object) {
            cat("An object of class \"",class(object),"\"\n",sep="")
            for (i in names(Meta(object))) {
              cat(i,"\n")
              print(Meta(object)[[i]])
            }
            print(dataTable(object))
          })

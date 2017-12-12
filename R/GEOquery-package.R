

#' Generic functions for GEOquery
#' 
#' The main documentation is in the Class documentation
#' 
#' @name GEOData-accessors
#' @aliases dataTable Accession Columns GPLList GSMList Meta Table
#' @author Sean Davis
#' @seealso \code{\link{GEOData-class}}
#' @keywords IO
NULL





#' Class "GDS"
#' 
#' A class describing a GEO GDS entity
#' 
#' 
#' @name GDS-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("GDS", ...)}
#' @author Sean Davis
#' @seealso \code{\link{GEOData-class}}
#' @keywords classes
NULL



#' Convert a GDS data structure to a BioConductor data structure
#' 
#' Functions to take a GDS data structure from getGEO and coerce it to limma
#' MALists or ExpressionSets.
#' 
#' This function just rearranges one data structure into another.  For GDS, it
#' also deals appropriately with making the "targets" list item for the limma
#' data structure and the phenoData slot of ExpressionSets.
#'
#' @name coercion
#' 
#' @aliases GDS2MA GDS2eSet
#' @param GDS The GDS datastructure returned by getGEO
#' @param do.log2 Boolean, should the data in the GDS be log2 transformed
#' before inserting into the new data structure
#' @param GPL Either a GPL data structure (from a call to getGEO) or NULL.  If
#' NULL, this will cause a call to getGEO to produce a GPL.  The gene
#' information from the GPL is then used to construct the \code{genes} slot of
#' the resulting limma \code{MAList} object or the \code{featureData} slot of
#' the \code{ExpressionSet} instance.
#' @param AnnotGPL In general, the annotation GPL files will be available for
#' GDS records, so the default is to use these files over the user-submitted
#' GPL files
#' @param getGPL A boolean defaulting to TRUE as to whether or not to download
#' and include GPL information when converting to ExpressionSet or MAList.  You
#' may want to set this to FALSE if you know that you are going to annotate
#' your featureData using Bioconductor tools rather than relying on information
#' provided through NCBI GEO.  Download times can also be greatly reduced by
#' specifying FALSE.
#' @return \item{GDS2MA}{A limma MAList} \item{GDS2eSet}{An ExpressionSet
#' object}
#' @author Sean Davis
#' @references See the limma and ExpressionSet help in the appropriate packages
#' @keywords IO
#' @examples
#' 
#' 
#' \dontrun{gds505 <- getGEO('GDS505')}
#' \dontrun{MA <- GDS2MA(gds505)}
#' \dontrun{eset <- GDS2eSet(gds505)}
#' 
#' 
NULL





#' Class "GEOData"
#' 
#' A virtual class for holding GEO samples, platforms, and datasets
#' 
#' 
#' @name GEOData-class
#' @aliases GEOData-class Accession,GEOData-method Columns,GEOData-method
#' Meta,GEOData-method Table,GEOData-method dataTable,GEOData-method
#' show,GEOData-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("GEOData", ...)}.
#' @author Sean Davis
#' @seealso \code{\link{GDS-class}}, \code{\link{GPL-class}},
#' \code{\link{GSM-class}}, \code{\link{GEODataTable-class}},
#' @keywords classes
NULL





#' Class "GEODataTable"
#' 
#' Contains the column descriptions and data for the datatable part of a GEO
#' object
#' 
#' 
#' @name GEODataTable-class
#' @aliases GEODataTable-class Accession,GEODataTable-method
#' Columns,GEODataTable-method Meta,GEODataTable-method
#' Table,GEODataTable-method dataTable,GEODataTable-method
#' show,GEODataTable-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("GEODataTable", ...)}.
#' @author Sean Davis
#' @keywords classes
NULL





#' Class "GPL"
#' 
#' Contains a full GEO Platform entity
#' 
#' @aliases GPL GPL,GDS-method
#' @name GPL-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("GPL", ...)}.
#' @author Sean Davis
#' @seealso \code{\link{GEOData-class}}
#' @keywords classes
NULL





#' Class "GSE"
#' 
#' Contains a GEO Series entity
#' 
#' 
#' @name GSE-class
#' @aliases GSE-class GPLList,GSE-method GSMList,GSE-method Meta,GSE-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("GSE", ...)}.
#' @author Sean Davis
#' @seealso \code{\link{GPL-class}},\code{\link{GSM-class}}
#' @keywords classes
NULL





#' Class "GSM"
#' 
#' A class containing a GEO Sample entity
#' 
#' 
#' @name GSM-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("GSM", ...)}.
#' @author Sean Davis
#' @seealso \code{\link{GEOData-class}}
#' @keywords classes
NULL




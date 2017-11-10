#' Get GSE data tables from GEO into R data structures.
#' 
#' In some cases, instead of individual sample records (GSM) containing
#' information regarding sample phenotypes, the GEO Series contains that
#' information in an attached data table.  And example is given by GSE3494
#' where there are two data tables with important information contained within
#' them.  Using getGEO with the standard parameters downloads the GSEMatrix
#' file which, unfortunately, does not contain the information in the data
#' tables.  This function simply downloads the ``header'' information from the
#' GSE record and parses out the data tables into R data.frames.
#' 
#' 
#' @param GSE The GSE identifier, such as ``GSE3494''.
#' @return A list of data.frames.
#' @author Sean Davis <sdavis2@@mail.nih.gov>
#' @seealso \code{\link{getGEO}}
#'
#' @importFrom xml2 xml_text xml_find_all read_xml
#' @importFrom readr read_tsv
#' 
#' @keywords IO
#' @examples
#' 
#' dfl = getGSEDataTables("GSE3494")
#' lapply(dfl,head)
#'
#' 
#' @export 
getGSEDataTables <- function(GSE) {
    url=sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&form=xml&view=full&acc=%s",GSE)
    doc1 = read_xml(url)
    dTableNodes = xml_find_all(doc1,"//d1:Data-Table")
    dTables = sapply(dTableNodes,function(x) {
        cnames=sapply(xml_find_all(x,"d1:Column/d1:Name"),xml_text)
        dTableText=xml_text(xml_find_all(x,"d1:Internal-Data")[[1]])
        browser()
        dTable = suppressWarnings(read_tsv(dTableText, col_names = FALSE))
        colnames(dTable)=cnames
        return(dTable)
    })
    return(dTables)
}

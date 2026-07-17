# Parse the <Data-Table> nodes of an acc.cgi form=xml Series document into a
# list of data.frames. Factored out of getGSEDataTables() so the parsing is
# unit-testable offline against a synthetic XML document (the only networked
# step is read_xml(url) in the caller). Using lapply (not sapply) guarantees a
# list of data.frames regardless of how many tables the Series carries -- sapply
# mis-simplifies a single table into a character matrix (#80).
.parseGSEDataTableNodes <- function(dTableNodes) {
    lapply(dTableNodes, function(x) {
        cnames = vapply(xml2::xml_find_all(x, "d1:Column/d1:Name"), xml2::xml_text, character(1))
        dTableText = xml2::xml_text(xml2::xml_find_all(x, "d1:Internal-Data")[[1]])
        dTable = suppressWarnings(readr::read_tsv(dTableText, col_names = FALSE))
        colnames(dTable) = cnames
        dTable
    })
}

#' Get GSE data tables from GEO into R data structures.
#'
#' In some cases, instead of individual sample records (GSM) containing
#' information regarding sample phenotypes, the GEO Series contains that
#' information in one or more data tables attached at the Series level. An
#' example is given by GSE3494, where there are two data tables with important
#' information contained within them. Series-level per-cell annotation tables
#' from single-cell studies (the `!series_table` blocks in SOFT format, e.g. the
#' "Listing of Individual Cells" table in GSE98638) are exposed here as well.
#' Using getGEO with the standard parameters downloads the GSEMatrix
#' file which, unfortunately, does not contain the information in the data
#' tables.  This function simply downloads the ``header'' information from the
#' GSE record and parses out the data tables into R data.frames.
#'
#'
#' @param GSE The GSE identifier, such as ``GSE3494''.
#' @return A list of data.frames, one per `<Data-Table>` in the Series record
#'   (a Series may carry zero, one, or several). Each data.frame's column names
#'   are taken from the table's column definitions.
#' @author Sean Davis <sdavis2@@mail.nih.gov>
#' @seealso \code{\link{getGEO}}
#'
#' @importFrom xml2 xml_text xml_find_all read_xml
#' @importFrom readr read_tsv
#'
#' @keywords IO
#' @examples
#' \dontrun{
#'
#' dfl = getGSEDataTables('GSE3494')
#' lapply(dfl,head)
#'
#'
#' }
#' @export
getGSEDataTables <- function(GSE) {
    url = sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&form=xml&view=full&acc=%s",
        GSE)
    doc1 = read_xml(url)
    dTableNodes = xml_find_all(doc1, "//d1:Data-Table")
    .parseGSEDataTableNodes(dTableNodes)
}

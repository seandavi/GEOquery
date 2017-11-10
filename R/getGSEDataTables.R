
# need xml2 readr
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

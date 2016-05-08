getGSEDataTables <- function(GSE) {
    url=sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&form=xml&view=full&acc=%s",GSE)
    ## Needed to switch to httr to support https, which is now what GEO is using
    txt = content(httr::GET(url),type='text')
    doc1 = xmlRoot(xmlParseDoc(txt))
    dTableNodes = getNodeSet(doc1,"//ns:Data-Table","ns")
    dTables = sapply(dTableNodes,function(x) {
	    cnames=sapply(getNodeSet(x,"ns:Column/ns:Name","ns"),xmlValue)
		dTableText=xmlValue(getNodeSet(x,"ns:Internal-Data","ns")[[1]])
		tc = textConnection(dTableText)
		dTable = read.delim(tc,sep="\t",header=FALSE)
		close(tc)
		colnames(dTable)=cnames
		return(dTable)
	})
    return(dTables)
}

getGSEDataTables <- function(GSE) {
    url=sprintf("http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?targ=self&form=xml&view=full&acc=%s",GSE)
    doc1 = xmlRoot(xmlParseDoc(url, options=HUGE, asText=FALSE))
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
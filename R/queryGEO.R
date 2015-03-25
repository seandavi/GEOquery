queryGEO = function(term,retmax=500,...) {
  tmp = entrez_summary(db='gds',id=entrez_search(db='gds',term,retmax=retmax,...)$ids)
  if(!is.list(tmp)) tmp = list(tmp)
  return(data.frame(do.call(rbind,lapply(tmp,function(y) {
    return(c(accession=y$Accession,
             taxon=y$taxon,
             entityType=y$entryType,
             pubDate = y$PDAT,
             nSamples=y$n_samples,
             title=y$title))
  }))))
}

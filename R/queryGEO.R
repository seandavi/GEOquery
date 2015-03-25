queryGEO = function(term,retmax=500,...) {
  tmp = list(entrez_summary(db='gds',id=entrez_search(db='gds',term,retmax=retmax,...)$ids))
  return(data.frame(do.call(rbind,lapply(tmp,function(y) {
    return(c(accession=y$Accession,
             taxon=y$taxon,
             entityType=y$entryType,
             pubDate = y$PDAT,
             nSamples=y$n_samples,
             title=y$title))
  }))))
}

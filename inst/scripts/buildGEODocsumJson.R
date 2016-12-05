# build geo docsum jsons
res = searchGEODocSums()

library(BiocParallel)
register(MulticoreParam())



docsumJSON = function(retstart,retmax) {
  require(jsonlite)
  dir.create(sprintf('~/docsums/docsums_%d',(retstart-1) %% 10000),showWarnings=FALSE,recursive = TRUE)
  f = gzfile(sprintf('~/docsums/docsums_%d/docsums_%d_%d.json.gz',(retstart-1) %% 10000,retstart,retstart+retmax),'wt')
  h = fetchDocSums(res,retstart=retstart,retmax=retmax)
  stream_out(h,f)
  message(retstart)
  close(f)
}

buildGEODocsums = function() {
  bplapply(seq(0,res$count,100),function(x) docsumJSON(x+1,100))
}

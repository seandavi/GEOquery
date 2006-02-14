"GDS2MA" <-
  function(GDS,do.log2=FALSE,GPL=NULL) {
    require(limma)
    if(is.null(GPL)) {
      GPL <- getGEO(Meta(GDS)$platform)
    }
    ord.table <- match(Table(GDS)[,1],Table(GPL)[,1])
                                        # exclude non-numeric columns
    inc.columns <- grep('GSM',colnames(Table(GDS)))
    if(do.log2) {
      M <- as.matrix(log2(Table(GDS)[ord.table,inc.columns]))
    } else {
      M <- as.matrix(Table(GDS)[ord.table,inc.columns])
    }
    MA <- new('MAList',list(M=M,
                            A=NULL,
                            targets=Columns(GDS),
                            genes=Table(GPL),
                            notes=Meta(GDS)
                            ))
    return(MA)
  }

"GDS2eSet" <-
  function(GDS,do.log2=FALSE) {
    require(Biobase)
                                        # exclude non-numeric columns
    inc.columns <- grep('GSM',colnames(Table(GDS)))
    if(do.log2) {
      expr <- as.matrix(log2(Table(GDS)[,inc.columns]))
    } else {
      expr <- as.matrix(Table(GDS)[,inc.columns])
    }
    rownames(expr) <- Table(GDS)$ID
    tmp <- Columns(GDS)
    rownames(tmp) <- as.character(tmp$sample)
    pheno <- new('phenoData',
                 pData=tmp,
                 varLabels=as.list(colnames(Columns(GDS))))
    eset <- new('exprSet',exprs=expr,phenoData=pheno)
    sampleNames
    geneNames(eset) <- Table(GDS)$ID_REF
    return(eset)
  }

"GDS2MA" <-
  function(GDS,do.log2=FALSE,GPL=NULL) {
    require(limma)
    if(is.null(GPL)) {
      GPL <- getGEO(Meta(GDS)$platform)
    }
    ord.table <- match(Table(GDS)[,1],Table(GPL)[,1])
                                        # exclude non-numeric columns
    inc.columns <- grep('GSM',colnames(Table(GDS)))
    mat <- suppressWarnings(as.matrix(apply(Table(GDS)[ord.table,inc.columns],2,function(x) {as.numeric(as.character(x))})))
    if(do.log2) {
      M <- log2(mat)
    } else {
      M <- mat
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
  function(GDS,do.log2=FALSE,GPL=NULL) {
    require(Biobase)
                                        # exclude non-numeric columns
    if(is.null(GPL)) {
      GPL <- getGEO(Meta(GDS)$platform)
    }
    ord.table <- match(Table(GDS)[,1],Table(GPL)[,1])
    inc.columns <- grep('GSM',colnames(Table(GDS)))
    mat <- suppressWarnings(as.matrix(apply(Table(GDS)[,inc.columns],2,
                                            function(x) {as.numeric(as.character(x))})))
    if(do.log2) {
      expr <- log2(mat)
    } else {
      expr <- mat
    }
    rownames(expr) <- Table(GDS)$ID_REF
    tmp <- Columns(GDS)
    rownames(tmp) <- as.character(tmp$sample)
    pheno <- new("AnnotatedDataFrame",data=tmp)
    mabstract=ifelse(is.null(Meta(GDS)$description),"",Meta(GDS)$description)
    mpubmedids=ifelse(is.null(Meta(GDS)$pubmed_id),"",Meta(GDS)$pubmed_id)
    mtitle=ifelse(is.null(Meta(GDS)$title),"",Meta(GDS)$title)
    featuredata <- new('AnnotatedDataFrame',data=Table(GPL),
                       varMetadata=data.frame(Column=Columns(GPL)[,1],
                         labelDescription=Columns(GPL)[,2]))
    eset <- new('ExpressionSet',exprs=expr,phenoData=pheno,
                featureData=featuredata,
                experimentData=new("MIAME",
                  abstract=mabstract,
                  title=mtitle,
                  pubMedIds=mpubmedids,
                  other=Meta(GDS)))
    return(eset)
  }

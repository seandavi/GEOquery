"GDS2MA" <-
  function(GDS,do.log2=FALSE,GPL=NULL,AnnotGPL=TRUE,getGPL=TRUE) {
    require(limma)
    if(!is.null(GPL) & getGPL) {
      GPL <- getGEO(Meta(GDS)$platform,AnnotGPL=AnnotGPL)
      ord.table <- match(Table(GDS)[,1],Table(GPL)[,1])
                                        # exclude non-numeric columns
    }
    inc.columns <- grep('GSM',colnames(Table(GDS)))
    mat <- suppressWarnings(as.matrix(apply(Table(GDS)[,inc.columns],2,function(x) {as.numeric(as.character(x))})))
    if(do.log2) {
      M <- log2(mat)
    } else {
      M <- mat
    }
    if(is.null(GPL)) {
        ## No GPL requested
        return(new('MAList',list(M=M,
                                 A=NULL,
                                 targets=Columns(GDS),
                                 notes=Meta(GDS)
                                 )))
    } else {
        ## include GPL
        return(new('MAList',list(M=M,
                                 A=NULL,
                                 targets=Columns(GDS),
                                 genes=Table(GPL)[ord.table,],
                                 notes=Meta(GDS)
                                 )))}

    return(MA)
  }

"GDS2eSet" <-
  function(GDS,do.log2=FALSE,GPL=NULL,AnnotGPL=TRUE,getGPL=TRUE) {
                                        # exclude non-numeric columns
    inc.columns <- grep('GSM',colnames(Table(GDS)))
    mat <- suppressWarnings(as.matrix(apply(Table(GDS)[,inc.columns],2,
                                            function(x) {as.numeric(as.character(x))})))
    if(do.log2) {
      expr <- log2(mat)
    } else {
      expr <- mat
    }
    rownames(expr) <- as.character(Table(GDS)$ID_REF)
    tmp <- Columns(GDS)
    rownames(tmp) <- as.character(tmp$sample)
    pheno <- new("AnnotatedDataFrame",data=tmp)
    mabstract=ifelse(is.null(Meta(GDS)$description),"",Meta(GDS)$description)
    mpubmedids=ifelse(is.null(Meta(GDS)$pubmed_id),"",Meta(GDS)$pubmed_id)
    mtitle=ifelse(is.null(Meta(GDS)$title),"",Meta(GDS)$title)
    featuredata = new("AnnotatedDataFrame",data=data.frame(row.names=rownames(expr)))
    if(getGPL | !is.null(GPL)) {
        if(is.null(GPL)) {
            GPL <- getGEO(Meta(GDS)$platform,AnnotGPL=AnnotGPL)
        }
        ord.table <- match(Table(GDS)[,1],Table(GPL)[,1])        
        dt <- Table(GPL)
        rownames(dt) <- as.character(dt$ID)
        vardt = data.frame(Column=Columns(GPL)[,1],
            labelDescription=Columns(GPL)[,2])
        ## GEO started using the same column names for
        ## both GO IDs and textual descriptions, so
        ## had to be made unique.
        vardt[,1] <- make.unique(as.character(vardt[,1]))
        rownames(vardt) = vardt[,1]
        colnames(dt) <- rownames(vardt)
        featuredata <- new('AnnotatedDataFrame',data=dt[ord.table,],
                           varMetadata=vardt)
    } 
    eset <- new('ExpressionSet',exprs=expr,phenoData=pheno,
                featureData=featuredata,
                experimentData=new("MIAME",
                  abstract=mabstract,
                  title=mtitle,
                  pubMedIds=mpubmedids,
                  other=Meta(GDS)))
    return(eset)
  }

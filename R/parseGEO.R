parseGEO <- function(fname,GSElimits) {
  con <- fileOpen(fname)
  first.entity <- findFirstEntity(con)
  close(con)
  ret <- switch(as.character(first.entity[1]),
                sample= {
                  parseGSM(fname)
                },
                series= parseGSE(fname,GSElimits),
                dataset= {
                  parseGDS(fname)
                },
                platform= {
                  parseGPL(fname)
                },
                "0" = {
                  parseGSEMatrix(fname)$eset
                },
                )
  return(ret)
}

parseGeoMeta <- function(txt) {
  ## leader <- strsplit(grep('!\\w*_',txt,perl=TRUE,value=TRUE)[1],'_')[[1]][1]
  ## pull out only lines that are in the header
  tmp <- txt[grep("!\\w*?_",txt)]
  tmp <- gsub("!\\w*?_",'',tmp)
  first.eq <- regexpr(' = ',tmp)
  tmp <- cbind(substring(tmp,first=1,last=first.eq-1),
               substring(tmp,first=first.eq+3))
  tmp <- tmp[tmp[,1]!="",]
  header <- split(tmp[,2],tmp[,1])
  return(header)
}

parseGDSSubsets <- function(txt) {
  ## takes GDS text as input
  ## returns a data frame suitable for inclusion as a "Column" slot
  ##   in a GDS GeoDataTable object
  numSubsets <- length(grep('^\\^subset',txt,ignore.case=TRUE))
  subset.lists <- list()
  if (numSubsets>0) {
    subset.types <-
      do.call('rbind',strsplit(txt[grep("!subset_type",txt)],' = '))[,2]
    subset.descriptions <-
      do.call('rbind',strsplit(txt[grep("!subset_description",txt)],' = '))[,2]
    subset.samples <-
      strsplit(do.call('rbind',strsplit(txt[grep("!subset_sample_id",txt)],' = '))[,2],',')
    for(i in 1:length(subset.types)) {
      for(j in 1:length(subset.samples[[i]])) {
        subset.lists[[subset.types[i]]][[subset.samples[[i]][j]]] <- subset.descriptions[i]
      }
    }
  }
  sample.descriptions <-
    do.call('rbind',strsplit(txt[grep('#GSM',txt)],' = '))
  sample.descriptions[,1] <- gsub('#','',sample.descriptions[,1])
  samp.table <- data.frame(sample=as.character(sample.descriptions[,1]))
  if(length(subset.lists)>0) {
    for(i in 1:length(subset.lists)) {
      samp.table[match(names(subset.lists[[i]]),samp.table[,1]),
                 names(subset.lists)[i]] <- as.factor(unlist(subset.lists[[i]]))
    }
    colnames(samp.table) <- c('sample',gsub(' ','.',names(subset.lists)))
  } else {
    colnames(samp.table) <- 'sample'
  }
  samp.table[match(sample.descriptions[,1],samp.table[,1]),'description'] <-
    sample.descriptions[,2]
  return(as.data.frame(samp.table))
}


splitOnFirst <- function(x,pattern) {
  patlen <- nchar(pattern)
  matches <- regexpr(pattern,x)
  leftside <- substr(x,start=1,stop=matches-1)
  rightside <- substr(x,start=matches+patlen,stop=10000000)
  return(data.frame(leftside,rightside))
}


parseGeoColumns <- function(txt) {
  cols <- as.data.frame(splitOnFirst(txt[grep('^#',txt,perl=TRUE)],' = '))
  cols[,1] <- sub('#','',as.character(cols[,1]))
  colnames(cols) <- c('Column','Description')
  return(cols)
}

### Parse a GSM; the limit is established by the
### parameter n, used by getGSE to limit the number
### of lines read to only the size of ONE GSM,
### since a single GSE contains many GSMs
.parseGSMWithLimits <- function(con,n=NULL) {
  txt <- vector('character')
  i <- 0
  hasDataTable=FALSE
  while(i <- i+1) {
    tmp <- try(readLines(con,1))
    if(inherits(tmp,"try-error") | length(tmp)==0) {
      hasDataTable=FALSE
      break
    }
    txt[i] <- tmp
    if(length(grep('!\\w+_table_begin',txt[i],perl=TRUE))>0) {
      hasDataTable=TRUE
      break
    }
  }
  cols <- parseGeoColumns(txt)
  meta <- parseGeoMeta(txt)
  geoDataTable <- new("GEODataTable",columns=data.frame(),table=data.frame())
  if(hasDataTable) {
    nLinesToRead <- NULL
    if(!is.null(n)) {
      nLinesToRead <- n-length(txt)
    }
    dat3 <- fastTabRead(con,n=nLinesToRead)
    geoDataTable <- new('GEODataTable',columns=cols,table=dat3[1:(nrow(dat3)-1),])
  } 
  gsm <- new('GSM',
             header=meta,
             dataTable = geoDataTable)
}

parseGSM <- function(fname) {
  con <- fileOpen(fname)
  ret <- .parseGSMWithLimits(con)
  close(con)
  return(ret)
}


### This function does a grep on a file
### by doing a readline in chunks of size
### chunksize.
### Return value is a data.frame with the line number
### of each match and the line itself.
filegrep <-
  function(con,regex,chunksize=10000) {
    i <- 0
    ret <- NULL
    while(TRUE) {
      lines <- readLines(con,n=chunksize)
      if(length(lines)==0) {
        break
      }
      foundLines <- grep(regex,lines)
      foundTypes <- lines[foundLines]
      if(length(foundLines)>0) {
        foundLines <- foundLines+i
        tmp <- data.frame(foundLines=foundLines,foundTypes=foundTypes)
        if(is.null(ret)) {
          ret <- tmp
        } else {
          ret <- rbind(ret,tmp)
        }
      }
      i <- i+length(lines)
    }
    return(ret) 
  }

parseGSE <- function(fname,GSElimits) {
  gsmlist <- list()
  gpllist <- list()
  GSMcount <- 0
  writeLines('Parsing....')
  con <- fileOpen(fname)
  lineCounts <- filegrep(con,"\\^(SAMPLE|PLATFORM)",chunksize=10000)
  cat(sprintf("Found %d entities...\n",nrow(lineCounts)))
  close(con)
  ## I close and reopen the file because on Windows, the seek
  ## function is pretty much guaranteed to NOT work
  con <- fileOpen(fname)
  ## This gets the header information for the GSE
  a <- readLines(con,lineCounts[1,1]-1)
  header=parseGeoMeta(a)
  ## parse the actual entities, now
  for(j in 1:nrow(lineCounts)) {
    tmp <- strsplit(as.character(lineCounts[j,2])," = ")[[1]]
    accession <- tmp[2]
    cat(sprintf("%s (%d of %d entities)\n",accession,j,nrow(lineCounts)))
    entityType <- tolower(sub("\\^","",tmp[1]))
    nLinesToRead <- lineCounts[j+1,1]-lineCounts[j,1]-1
    if(j==nrow(lineCounts)) {
      nLinesToRead <- NULL
    }
    if(entityType=="sample") {
      gsmlist[[accession]] <- .parseGSMWithLimits(con,n=nLinesToRead)
    }
    if(entityType=="platform") {
      gpllist[[accession]] <- .parseGPLWithLimits(con,n=nLinesToRead)
    }
  }
  close(con)
  return(new("GSE",
             header= header,
             gsms  = gsmlist,
             gpls  = gpllist))
}


findFirstEntity <- function(con) {
  while(TRUE) {
    line <- readLines(con,1)
    if(length(line)==0) return(0)
    entity.line <- grep('^\\^(DATASET|SAMPLE|SERIES|PLATFORM|ANNOTATION)',
                        line,ignore.case=TRUE,value=TRUE,perl=TRUE)
    entity.line <- gsub('annotation','platform',entity.line,ignore.case=TRUE)
    if(length(entity.line)>0) {
      ret <- c(tolower(sub('\\^','',strsplit(entity.line,' = ')[[1]][1])),
               strsplit(entity.line,' = ')[[1]][2])
      return(ret)
    }
  }
}

fastTabRead <- function(con,sep="\t",header=TRUE,sampleRows=100,
                        colClasses=NULL,n=NULL,...) {
### Need to read tables quickly, so guess the colclasses on the
### fly.  This is a bit dangerous since the first rows might not
### be representative of the entire set, but it is SO MUCH FASTER
### than the alternative, I have to do it.
  dat3 <- data.frame(NULL)
  numberOfLines <- -1
  if(!is.null(n)) {
    numberOfLines <- n-sampleRows
  }
  if(is.null(colClasses)) {
    if(!is.null(n)) {
      sampleRows <- min(sampleRows,n)
    }
    dat1 <- read.delim(con,sep=sep,header=TRUE,nrows=sampleRows,quote="",comment.char="",na.strings=c('NA','null','NULL'),...)
    colclasses <- apply(dat1,2,class)
    colclasses[1] <- "factor"
    dat2 <- read.delim(con,sep=sep,colClasses=colclasses,
                       header=FALSE,quote="",comment.char="",
                       na.strings=c('NA','null','NULL'),
                       nrows=numberOfLines,...)
    colnames(dat2) <- colnames(dat1)
    dat3 <- rbind(dat1,dat2)
  } else {
    dat3 <- read.delim(con,sep=sep,colClasses=colClasses,
                       header=header,quote="",comment.char="",
                       na.strings=c('NA','null','NULL'),nrows=numberOfLines,...)
  }
  return(dat3)
}

parseGDS <- function(fname) {
  con <- fileOpen(fname)
  txt <- vector('character')
  i <- 0
  while(i <- i+1) {
    txt[i] <- readLines(con,1)
    if(length(grep('!\\w+_table_begin',txt[i],perl=TRUE))>0) break
  }
  cols <- parseGDSSubsets(txt)
  meta <- parseGeoMeta(txt)
  dat3 <- fastTabRead(con)
  close(con)
  geoDataTable <- new('GEODataTable',columns=cols,table=dat3[1:(nrow(dat3)-1),])
  gds <- new('GDS',
             header=meta,
             dataTable = geoDataTable)
  return(gds)
}


.parseGPLWithLimits <- function(con,n=NULL) {
  txt <- vector('character')
  i <- 0
  hasDataTable=FALSE
  while(i <- i+1) {
    tmp <- try(readLines(con,1))
    if(inherits(tmp,"try-error") | length(tmp)==0 ) {
      hasDataTable=FALSE
      break
    }
    if(!is.null(n)) {
      if(i==n) {
        hasDataTable=FALSE
        break
      }
    }
    txt[i] <- tmp
    if(length(grep('!\\w+_table_begin',txt[i],perl=TRUE))>0) {
      hasDataTable=TRUE
      break
    }
  }
  cols <- parseGeoColumns(txt)
  meta <- parseGeoMeta(txt)
  geoDataTable <- new("GEODataTable",columns=data.frame(),table=data.frame())
  if(hasDataTable) {
    nLinesToRead <- NULL
    if(!is.null(n)) {
      nLinesToRead <- n-length(txt)
    }
    dat3 <- fastTabRead(con,n=nLinesToRead)
    geoDataTable <- new('GEODataTable',columns=cols,table=dat3[1:(nrow(dat3)-1),])
  } 
  gpl <- new('GPL',
             header=meta,
             dataTable = geoDataTable)
}

parseGPL <- function(fname) {
  con <- fileOpen(fname)
  ret <- .parseGPLWithLimits(con)
  close(con)
  return(ret)
}

txtGrab <- function(regex,x) {
  x <- as.character(x)
  a <- regexpr(regex,x,perl=TRUE)
  return(substr(x,a,a+attr(a,'match.length')-1))
}

### Function wrapper to get and parse ALL
### the GSEMatrix files associated with a GSE
### into a list of ExpressionSets
getAndParseGSEMatrices <- function(GEO,destdir) {
  GEO <- toupper(GEO)
  ## This stuff functions to get the listing of available files
  ## for a given GSE given that there may be many GSEMatrix
  ## files for a given GSE.
  a <- getURL(sprintf('ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/%s/',GEO))
  tmpcon <- textConnection(a,'r')
  b <- read.table(tmpcon)
  close(tmpcon)
  b <- as.character(b[,ncol(b)])
  message(sprintf('Found %d file(s)',length(b)))
  ret <- list()
  ## Loop over the files, returning a list, one element
  ## for each file
  for(i in 1:length(b)) {
    message(b[i])
    destfile=file.path(destdir,b[i])
    if(file.exists(destfile)) {
      message(sprintf('Using locally cached version: %s',destfile))
    } else {
      download.file(sprintf('ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/%s/%s',
                            GEO,b[i]),destfile=destfile,mode='wb')
    }
    ret[[b[i]]] <- parseGSEMatrix(destfile)$eset
  }
  return(ret)
}


### Function to parse a single GSEMatrix
### file into an ExpressionSet
parseGSEMatrix <- function(fname) {
  require(Biobase)
  dat <- readLines(fname)
  ## get the number of !Series and !Sample lines
  nseries <- sum(grepl("^!Series_", dat))
  nsamples <- sum(grepl("^!Sample_", dat))
  con <- fileOpen(fname)
  ## Read the !Series_ and !Sample_ lines
  header <- read.table(con,sep="\t",header=FALSE,nrows=nseries)
  tmpdat <- read.table(con,sep="\t",header=FALSE,nrows=nsamples)
  sampledat <- data.frame(t(tmpdat[,-1]))
  colnames(sampledat) <- make.unique(sub('!Sample_','',as.character(tmpdat[,1])))
  readLines(con,1)
  colClasses <- c('character',rep('numeric',nrow(sampledat)))
  datamat <- as.matrix(read.delim(con,sep="\t",header=TRUE,row.names=1,
                                  colClasses=colClasses,
                                  na.strings=c('NA','null','NULL'),
                                  comment.char=""))
  close(con)
  ## All the series matrix files are assumed to end with
  ## the line "!series_matrix_table_end", so we remove
  ## that line from the datamatrix (it has been read)
  datamat <- datamat[1:(nrow(datamat)-1),]
  rownames(sampledat) <- colnames(datamat)
  GPL=as.character(sampledat[1,grep('platform_id',colnames(sampledat),ignore.case=TRUE)])
  gpl <- getGEO(GPL)
  vmd <- Columns(gpl)
  rownames(vmd) <- colnames(Table(gpl))
  dat <- Table(gpl)
  rownames(dat) <- as.character(dat$ID)
  dat <- dat[match(rownames(datamat),rownames(dat)),]
  fd <- new('AnnotatedDataFrame',data=dat,varMetadata=vmd)
  if(is.null(nrow(datamat))) {
    datamat=matrix(nrow=0,ncol=nrow(sampledat))
  }
  eset <- new('ExpressionSet',
              phenoData=as(sampledat,'AnnotatedDataFrame'),
              annotation=GPL,
              featureData=fd,
              exprs=datamat)
  return(list(GPL=as.character(sampledat[1,grep('platform_id',colnames(sampledat),ignore.case=TRUE)]),eset=eset))
}

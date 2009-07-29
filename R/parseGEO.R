parseGEO <- function(con,GSElimits) {
  first.entity <- findFirstEntity(con)
  ret <- switch(as.character(first.entity[1]),
                sample= {
                  parseGSM(con)
                },
                series= parseGSE(con,GSElimits),
                dataset= {
                  parseGDS(con)
                },
                platform= {
                  parseGPL(con)
                },
                )
  return(ret)
}

parseGeoMeta <- function(txt) {
		
	leader <- strsplit(grep('!\\w*_',txt,perl=TRUE,value=TRUE)[1],'_')[[1]][1]
     # pull out only lines that are in the header
  	tmp <- txt[grep(leader,txt)]
  	tmp <- gsub(paste(leader,'_',sep=""),'',tmp)
  	first.eq <- regexpr(' = ',tmp)
  	tmp <- cbind(substring(tmp,first=1,last=first.eq-1),substring			    (tmp,first=first.eq+3))
     # remove blank lines
     #tmp <- tmp[-grep('^\\s?$',tmp[,2],perl=TRUE),]
	header <- split(tmp[,2],tmp[,1])
	return(header)
}

parseGDSSubsets <- function(txt) {
  # takes GDS text as input
  # returns a data frame suitable for inclusion as a "Column" slot
  #   in a GDS GeoDataTable object
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

.parseGSMWithLimits <- function(con,n=NULL) {
  txt <- vector('character')
  i <- 0
  while(i <- i+1) {
    txt[i] <- readLines(con,1)
    if(length(grep('!\\w+_table_begin',txt[i],perl=TRUE))>0) break
  }
  cols <- parseGeoColumns(txt)
  meta <- parseGeoMeta(txt)
  nLinesToRead <- NULL
  if(!is.null(n)) {
    nLinesToRead <- n-length(txt)
  }
  dat3 <- fastTabRead(con,n=nLinesToRead)
  geoDataTable <- new('GEODataTable',columns=cols,table=dat3[1:(nrow(dat3)-1),])
  gsm <- new('GSM',
             header=meta,
             dataTable = geoDataTable)
}

parseGSM <- function(con) {
  return(.parseGSMWithLimits(con))
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

parseGSE <- function(con,GSElimits) {
  gsmlist <- list()
  gpllist <- list()
  GSMcount <- 0
  writeLines('Parsing....')
  lineCounts <- filegrep(con,"\\^(SAMPLE|PLATFORM)",chunksize=10000)
  cat(sprintf("Found %d entities...\n",nrow(lineCounts)))
  seek(con,0)
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

parseGDS <- function(con) {
    txt <- vector('character')
    i <- 0
    while(i <- i+1) {
      txt[i] <- readLines(con,1)
      if(length(grep('!\\w+_table_begin',txt[i],perl=TRUE))>0) break
    }
    cols <- parseGDSSubsets(txt)
    meta <- parseGeoMeta(txt)
    dat3 <- fastTabRead(con)
    geoDataTable <- new('GEODataTable',columns=cols,table=dat3[1:(nrow(dat3)-1),])
    gds <- new('GDS',
               header=meta,
               dataTable = geoDataTable)
  }


.parseGPLWithLimits <- function(con,n=NULL) {
  txt <- vector('character')
  i <- 0
  while(i <- i+1) {
    txt[i] <- readLines(con,1)
    if(length(grep('!\\w+_table_begin',txt[i],perl=TRUE))>0) break
  }
  cols <- parseGeoColumns(txt)
  meta <- parseGeoMeta(txt)
  nLinesToRead <- NULL
  if(!is.null(n)) {
    nLinesToRead <- n-length(txt)
  }
  dat3 <- fastTabRead(con,n=nLinesToRead)
  geoDataTable <- new('GEODataTable',columns=cols,table=dat3[1:(nrow(dat3)-1),])
  gpl <- new('GPL',
             header=meta,
             dataTable = geoDataTable)
}

parseGPL <- function(con) {
  return(.parseGPLWithLimits(con))
}

txtGrab <- function(regex,x) {
  x <- as.character(x)
  a <- regexpr(regex,x,perl=TRUE)
  return(substr(x,a,a+attr(a,'match.length')-1))
}

getAndParseGSEMatrices <- function(GEO) {
  require(RCurl)
  GEO <- toupper(GEO)
  a <- getURL(sprintf('ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/%s/',GEO))
  tmpcon <- textConnection(a,'r')
  b <- read.table(tmpcon)
  close(tmpcon)
  b <- as.character(b[,ncol(b)])
  writeLines(sprintf('Found %d file(s)',length(b)))
  ret <- list()
  for(i in 1:length(b)) {
    writeLines(b[i])
    tmp <- tempdir()
    download.file(sprintf('ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/%s/%s',GEO,b[i]),destfile=file.path(tmp,b[i]),mode='wb')
    tmpfile <- tempfile()
    gunzip(file.path(tmp,b[i]),tmpfile)
    con <- file(tmpfile,'r')
    ret[[b[i]]] <- parseGSEMatrix(con)$eset
    close(con)
  }
  return(ret)
}
    
  

parseGSEMatrix <- function(con) {
  require(Biobase)
  i <- 0
  while(i <- i+1) {
    a <- readLines(con,1)
    if(length(grep('^!Series_',a,ignore.case=TRUE))==0) {
      break
    }
  }
  seek(con,where=0)
  header <- read.table(con,sep="\t",header=FALSE,nrows=i-1)
  i <- 0
  seekloc <- 0
  while(1) {
    a <- readLines(con,1)
    if(length(grep('^!Sample_',a,ignore.case=TRUE))==0) {
      if(i==0) {
        seekloc <- seek(con)
        next
      } else {
        break
      }
    }
    i <- i+1
  }
  seek(con,seekloc)
  tmpdat <- read.table(con,sep="\t",header=FALSE,nrows=i)
  sampledat <- data.frame(t(tmpdat[,-1]))
  colnames(sampledat) <- make.unique(sub('!Sample_','',as.character(tmpdat[,1])))
  readLines(con,1)
  colClasses <- c('character',rep('numeric',nrow(sampledat)))
  datamat <- as.matrix(read.delim(con,sep="\t",header=TRUE,row.names=1,
                                  colClasses=colClasses,
                                  na.strings=c('NA','null','NULL'),
                                  comment.char=""))
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
  eset <- new('ExpressionSet',
              phenoData=as(sampledat,'AnnotatedDataFrame'),
              annotation=GPL,
              featureData=fd,
              exprs=datamat)
  return(list(GPL=as.character(sampledat[1,grep('platform_id',colnames(sampledat),ignore.case=TRUE)]),eset=eset))
}

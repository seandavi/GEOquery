parseGeoData <- function(txt) {
	# Get the stuff between table_begin and table_end
	# Arguments:
	#    txt: text associated with a single entity
	# Returns:
	#    data.frame containing just the table and associated
	#    column names
		
	# Get start and end of data table
	# Thanks to GEO folks for adding table_begin and end
	tbl.start <- grep('!\\w+_table_begin',txt,perl=TRUE)+1
	tbl.end <- grep('!\\w+_table_end',txt,perl=TRUE)-1
        if(length(tbl.end)==0) {
          tbl.end <- length(txt)
        }
#       This is slower, but correctly deals with column types
        txtcon <- textConnection(txt[tbl.start:tbl.end])
        ret <- read.delim(txtcon,header=TRUE,sep="\t",na.strings="NULL")
        close(txtcon)
        return(ret)
	# Grab first line after table_begin for column names
#	tbl.colnames <- strsplit(txt[tbl.start],"\t")[[1]]
#	tbl.tmp <- do.call('rbind',strsplit(txt[(tbl.start+1):tbl.end],"\t"))
#	colnames(tbl.tmp) <- tbl.colnames
#	return(data.frame(tbl.tmp))
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

parseGeoDataTable <- function(txt) {
	cols <- parseGeoColumns(txt)
	tab <- parseGeoData(txt)
	return(new('GEODataTable',columns=cols,table=tab))
}

parseGSM <- function(txt) {
  geoDataTable <- parseGeoDataTable(txt)
  meta <- parseGeoMeta(txt)
  gsm <- new('GSM',
             header=meta,
             dataTable = geoDataTable)
  return(gsm)
}

parseGPL <- function(txt) {
	geoDataTable <- parseGeoDataTable(txt)
	meta <- parseGeoMeta(txt)
	gsm <- new('GPL',
			  header=meta,
			  dataTable = geoDataTable)
		
	return(gsm)
}

parseGSE <- function(con,GSElimits) {
  gsmlist <- list()
  gpllist <- list()
  GSMcount <- 0
  writeLines('Parsing....')
  ## This gets the header information for the GSE
  lines <- 1
  a <- vector()
  finished <- FALSE
  nextEntity <- ""
  while(!finished) {
    line <- readLines(con,1)
    if(length(line)==0) finished <- TRUE
    a[lines] <- line
    lines <- lines+1
    b <- grep('^\\^(SAMPLE|PLATFORM)',line,value=TRUE,perl=TRUE)
    if(length(b)>0) {
      nextEntity <- b
      writeLines(b)
      finished <- TRUE
      lines <- 1
      header=parseGeoMeta(a)
    }
  }
  finished <- FALSE
  while(!finished) {
    line <- readLines(con,1)
    if(length(line)==0) {
      finished <- TRUE
    } else {
      a[lines] <- line
      lines <- lines+1
      b <- grep('^\\^(SAMPLE|PLATFORM)',line,value=TRUE,perl=TRUE)
    }
    if(length(b)>0 | finished) {
      lines <- 1
      #new SAMPLE
      if(length(grep('SAMPLE',nextEntity))>0) {
        accession <- strsplit(nextEntity,' = ')[[1]][2]
        GSMcount <- GSMcount+1
        # Look to see if limits should be used, otherwise, proceed
        if(is.null(GSElimits)) {
          gsmlist[[accession]] <- parseGSM(a)
          writeLines(b)
        } else {
          if((GSMcount>=GSElimits[1]) &
             (GSMcount<=GSElimits[2])) {
            gsmlist[[accession]] <- parseGSM(a)
            writeLines(b)
          } else {
            cat('Skipping sample',GSMcount,': Accession',accession,'at user request\n')
          }
        }
      }
      if(length(grep('PLATFORM',nextEntity))>0) {
        accession <- strsplit(nextEntity,' = ')[[1]][2]
        gpllist[[accession]] <- parseGPL(a)
      }
      if(!is.null(GSElimits)) {
        if(GSMcount+1>GSElimits[2]) {
                                        # end if beyond GSElimits[2]
          cat('Stopping here at user request\n')
          finished <- TRUE
        }
      }
      if(!finished) {
        nextEntity <- b
      }
      a <- vector()
    }
  }
  gse <- new("GSE",
             header=header,
             gsms = gsmlist,
             gpls = gpllist
             )
  return(gse)
}


parseGDS <- function(txt) {
  writeLines("parsing geodata")
  tab <- parseGeoData(txt)
  writeLines("parsing subsets")
  cols <- parseGDSSubsets(txt)
  
  writeLines("ready to return")
  return(new('GDS',
             header=parseGeoMeta(txt),
             dataTable=new('GEODataTable',
               table=tab,
               columns=cols
               )
             ))
	}

findFirstEntity <- function(con) {
  while(TRUE) {
    line <- readLines(con,1)
    if(length(line)==0) return(0)
    entity.line <- grep('^\\^(DATASET|SAMPLE|SERIES|PLATFORM)',
                        line,ignore.case=TRUE,value=TRUE,perl=TRUE)
    if(length(entity.line)>0) {
      ret <- c(tolower(sub('\\^','',strsplit(entity.line,' = ')[[1]][1])),
               strsplit(entity.line,' = ')[[1]][2])
      return(ret)
    }
  }
}

parseGPL2 <- function(con) {
  txt <- vector('character')
  i <- 0
  while(i <- i+1) {
    txt[i] <- readLines(con,1)
    if(length(grep('!\\w+_table_begin',txt[i],perl=TRUE))>0) break
  }
  cols <- parseGeoColumns(txt)
  meta <- parseGeoMeta(txt)
  dat1 <- read.delim(con,sep="\t",header=TRUE,nrows=10)
  colclasses <- apply(dat1,2,class)
  dat2 <- read.delim(con,sep="\t",colClasses=colclasses,header=FALSE)
  colnames(dat2) <- colnames(dat1)
  dat3 <- rbind(dat1,dat2)
  geoDataTable <- new('GEODataTable',columns=cols,table=dat3[1:(nrow(dat3)-1),])
  gsm <- new('GPL',
             header=meta,
             dataTable = geoDataTable)
}

parseGEO <- function(con,GSElimits) {
  first.entity <- findFirstEntity(con)
  ret <- switch(as.character(first.entity[1]),
                sample= {
                  txt <- readLines(con)
                  parseGSM(txt)
                },
                series= parseGSE(con,GSElimits),
                dataset= {
                  txt <- readLines(con)
                  parseGDS(txt)
                },
                platform= {
                  parseGPL2(con)
                },
                )
  return(ret)
}

txtGrab <- function(regex,x) {
  x <- as.character(x)
  a <- regexpr(regex,x,perl=TRUE)
  return(substr(x,a,a+attr(a,'match.length')-1))
}

getAndParseGSEMatrices <- function(GEO) {
  require(RCurl)
  require(R.utils)
  GEO <- toupper(GEO)
  a <- getURL(sprintf('ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/%s/',GEO))
  b <- read.table(textConnection(a,'r'))
  b <- as.character(b[,ncol(b)])
  writeLines(sprintf('Found %d file(s)',length(b)))
  ret <- list()
  for(i in 1:length(b)) {
    writeLines(b[i])
    tmp <- tempdir()
    download.file(sprintf('ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/%s/%s',GEO,b[i]),destfile=file.path(tmp,b[i]))
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
                                  colClasses=colClasses,comment.char=""))
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

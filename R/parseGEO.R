.na_strings = c('NA','null','NULL','Null')



#' Parse GEO text
#' 
#' Workhorse GEO parsers.
#' 
#' These are probably not useful to the end-user.  Use getGEO to access these
#' functions.  parseGEO simply delegates to the appropriate specific parser.
#' There should be no reason to use the parseGPL, parseGDS, parseGSE, or
#' parseGSM functions directly.
#' 
#' @aliases parseGEO parseGPL parseGSE parseGDS parseGSM
#' @param fname The filename of a SOFT format file.  If the filename ends in
#' .gz, a gzfile() connection is used to read the file directly.
#' @param GSElimits Used to limit the number of GSMs parsed into the GSE
#' object; useful for memory management for large GSEs.
#' @param destdir The destination directory into which files will be saved (to
#' be used for caching)
#' @param AnnotGPL Fetch the annotation GPL if available
#' @param getGPL Fetch the GPL associated with a GSEMatrix entity (should
#' remain TRUE for all normal use cases)
#' @return parseGEO returns an object of the associated type.  For example, if
#' it is passed the text from a GDS entry, a GDS object is returned.
#' @author Sean Davis
#' @seealso \code{\link{getGEO}}
#' @keywords IO
#' 
#' @export
parseGEO <- function(fname,GSElimits,destdir=tempdir(),AnnotGPL=FALSE,getGPL=TRUE) {
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
                      parseGSEMatrix(fname,destdir=destdir,AnnotGPL=AnnotGPL,getGPL=getGPL)$eset
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
#' @importFrom readr read_lines
.parseGSMWithLimits <- function(con,n=NULL) {
    txt <- vector('character')
    i <- 0
    hasDataTable=FALSE
    while(i <- i+1) {
        tmp <- try(read_lines(con,1))
        if(inherits(tmp,"try-error") | length(tmp)==0) {
            i <- i-1  # counted an extra line
            hasDataTable=FALSE
            break
        }
        txt[i] <- tmp
        if(i==length(txt)) txt <- c(txt, character(i))  # double vector size
        if(length(grep('!\\w+_table_begin',txt[i],perl=TRUE))>0) {
            hasDataTable=TRUE
            break
        }
    }
    txt <- txt[1:i]
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

#' @importFrom readr read_lines
parseGSE <- function(fname,GSElimits=NULL) {
    gsmlist <- list()
    gpllist <- list()
    GSMcount <- 0
    message('Reading file....')
    lines = read_lines(fname)
    message('Parsing....')
    lineCounts <- grep("\\^(SAMPLE|PLATFORM)",lines)
    message(sprintf("Found %d entities...",length(lineCounts)))
    ## I close and reopen the file because on Windows, the seek
    ## function is pretty much guaranteed to NOT work
    ## This gets the header information for the GSE
    a <- lines[1:(lineCounts[1]-1)]
    header=parseGeoMeta(a)
    lineCounts = c(lineCounts,length(lines))
    ## parse the actual entities, now
    for(j in 1:(length(lineCounts)-1)) {
        tmp <- strsplit(lines[lineCounts[j]]," = ")[[1]]
        accession <- tmp[2]
        message(sprintf("%s (%d of %d entities)",accession,j,length(lineCounts)))
        entityType <- tolower(sub("\\^","",tmp[1]))
        if(entityType=="sample") {
            gsmlist[[accession]] <- .parseGSMTxt(lines[lineCounts[j]:(lineCounts[j+1]-1)])
        }
        if(entityType=="platform") {
            gpllist[[accession]] <- .parseGPLTxt(lines[lineCounts[j]:(lineCounts[j+1]-1)])
        }
    }
    return(new("GSE",
               header= header,
               gsms  = gsmlist,
               gpls  = gpllist))
}


findFirstEntity <- function(con) {
    while(TRUE) {
        line <- suppressWarnings(readLines(con,100))
        entity.line <- grep('^\\^(DATASET|SAMPLE|SERIES|PLATFORM|ANNOTATION)',
                            line,ignore.case=TRUE,value=TRUE,perl=TRUE)
        entity.line <- gsub('annotation','platform',entity.line,ignore.case=TRUE)
        if(length(entity.line)>0) {
            ret <- c(tolower(sub('\\^','',strsplit(entity.line,' = ')[[1]][1])),
                     strsplit(entity.line,' = ')[[1]][2])
            return(ret)
        }
        # catch GSE SeriesMatrix files
        checkline = grep('^!Series_title',line,ignore.case=TRUE,value=TRUE,perl=TRUE)
        # Messy, but GSEMatrix files use tab-separation rather
        # than '=' for separation, so if the ' = ' is not present,
        # we presume to have a GSEMatrix file (return 0)
        if(length(checkline) > 0) {
            if(!grepl(' = ',checkline)) return(0)
        }
    }
    return(0)
}

fastTabRead <- function(con,sep="\t",header=TRUE,sampleRows=100,
                        colClasses=NULL,n=NULL,quote='"',...) {
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
        dat1 <- read.table(con,sep=sep,header=header,nrows=sampleRows,fill=TRUE,check.names=FALSE,
                           quote=quote,comment.char="",na.strings=c('NA','null','NULL','Null'),...)
        colclasses <- apply(dat1,2,class)
        colclasses[1] <- "factor"
        dat2 <- try(read.delim(con,sep=sep,colClasses=colclasses,
                               header=FALSE,quote=quote,comment.char="",
                               na.strings=c('NA','null','NULL','Null'),
                               nrows=numberOfLines,...),silent=TRUE)
        if(inherits(dat2,'try-error')) {
            dat3=dat1
        } else {
            colnames(dat2) <- colnames(dat1)
            dat3 <- rbind(dat1,dat2)
        }
    } else {
        dat3 <- read.delim(con,sep=sep,colClasses=colClasses,
                           header=header,quote=quote,comment.char="",
                           na.strings=c('NA','null','NULL',"Null"),nrows=numberOfLines,...)
    }
    return(dat3)
}



#' parse a GEO dataset (GDS)
#' 
#' GEO datasets, or GDS, used to be produced by NCBI GEO as "curated" versions
#' of submitted data. NCBI GEO no longer produces these objects, but there are
#' still thousands available that, in some cases, have much nicer sample
#' annotation and more standardized annotation for the associated GPLs.
#'
#' @importFrom readr read_lines read_tsv
#' 
#' @param fname the filename of the SOFT format file. May be gzipped.
#' @keywords internal
parseGDS <- function(fname) {
    txt = read_lines(fname)
    tbl_begin = grep('!\\w+_table_begin',txt,perl=TRUE)
    if(length(tbl_begin>0)) {
        txt = txt[1:tbl_begin[1]]
        dat3 <- read_tsv(fname, comment='!dataset_table_end', skip = tbl_begin[1],
                         guess_max = 10000, na = .na_strings)
    }
    cols <- parseGDSSubsets(txt)
    meta <- parseGeoMeta(txt)
    geoDataTable <- new('GEODataTable',columns=cols,table=as.data.frame(dat3))
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
        if(inherits(tmp,"try-error")
           || length(tmp)==0
           || (!is.null(n) && i==n)) {
            i <- i-1  # counted an extra line
            hasDataTable=FALSE
            break
        }
        txt[i] <- tmp
        if(i==length(txt)) txt <- c(txt, character(i))  # double vector size
        if(length(grep('!\\w+_table_begin',txt[i],perl=TRUE))>0) {
            hasDataTable=TRUE
            break
        }
    }
    txt <- txt[1:i]
    cols <- parseGeoColumns(txt)
    meta <- parseGeoMeta(txt)
    geoDataTable <- new("GEODataTable",columns=data.frame(),table=data.frame())
    if(hasDataTable) {
        nLinesToRead <- NULL
        if(!is.null(n)) {
            nLinesToRead <- n-length(txt)
        }
        dat3 <- fastTabRead(con,n=nLinesToRead,quote='')
        geoDataTable <- new('GEODataTable',columns=cols,table=dat3[1:(nrow(dat3)-1),])
    } 
    gpl <- new('GPL',
               header=meta,
               dataTable = geoDataTable)
}


#' @importFrom readr read_tsv
.parseGSMTxt <- function(txt) {
    tbl_begin = grep('!\\w+_table_begin',txt,perl=TRUE)
    if(length(tbl_begin>0)) {
        tbltxt = txt[(tbl_begin[1]+1):length(txt)]
        txt = txt[1:tbl_begin[1]]
        dat3 <- read_tsv(paste(tbltxt,collapse="\n"), comment='!sample_table_end',
                         guess_max = 10000, na = .na_strings)
    } else {
        # empty data table
        dat3 = data.frame()
    }

    cols <- parseGeoColumns(txt)
    meta <- parseGeoMeta(txt)
    geoDataTable <- new('GEODataTable',columns=cols,table=as.data.frame(dat3))
    geo <- new('GSM',
               header=meta,
               dataTable = geoDataTable)
    return(geo)
}
    

#' @importFrom readr read_lines
parseGSM <- function(fname) {
    txt = read_lines(fname)
    # read_lines reads separate blank lines
    # on windows, so remove them before 
    # proceeding. NOT doing so results in 
    # the Table header being repeated as the
    # first line of the Table and test failures
    # galore.
    txt = txt[txt != '']
    return(.parseGSMTxt(txt))
}

### In memory cache for GPL objects parsed from locally cached versions of GPL SOFT files.
### It is disabled by default with options('GEOquery.inmemory.gpl'=FALSE).
GPLcache <- new.env(parent=emptyenv())


#' @importFrom readr read_tsv
.parseGPLTxt <- function(txt) {
    tbl_begin = grep('!\\w+_table_begin',txt,perl=TRUE)
    
    if(length(tbl_begin>0)) {
        tbltxt = txt[(tbl_begin[1]+1):length(txt)]
        txt = txt[1:tbl_begin[1]]
        dat3 <- suppressMessages(read_tsv(paste(tbltxt,collapse="\n"), comment='!platform_table_end',
                         guess_max = 10000, na = .na_strings))
    } else {
        # empty data table
        dat3 = data.frame()
    }
    cols <- parseGeoColumns(txt)
    meta <- parseGeoMeta(txt)
    geoDataTable <- new('GEODataTable',columns=cols,table=as.data.frame(dat3))
    geo <- new('GPL',
               header=meta,
               dataTable = geoDataTable)
    return(geo)
}


#' @importFrom readr read_lines
parseGPL <- function(fname) {
    if(getOption('GEOquery.inmemory.gpl')) {
        info <- file.info(fname,extra_cols=FALSE)
        cache <- get0(fname,envir=GPLcache,inherits=FALSE)
        ## Check if the locally cached version wasn't modified.
        if(!is.null(cache) && cache$info$mtime==info$mtime) {
            message("Using GPL object found in memory from locally cached version.")
            return(cache$gpl)
        }
    }
    txt = read_lines(fname)
    # read_lines reads separate blank lines
    # on windows, so remove them before 
    # proceeding. NOT doing so results in 
    # the Table header being repeated as the
    # first line of the Table and test failures
    # galore.
    txt = txt[txt != '']
    return(.parseGPLTxt(txt))
}
    

txtGrab <- function(regex,x) {
    x <- as.character(x)
    a <- regexpr(regex,x,perl=TRUE)
    return(substr(x,a,a+attr(a,'match.length')-1))
}

### Function wrapper to get and parse ALL
### the GSEMatrix files associated with a GSE
### into a list of ExpressionSets
getAndParseGSEMatrices <- function(GEO,destdir,AnnotGPL,getGPL=TRUE,parseCharacteristics=TRUE) {
    GEO <- toupper(GEO)
    ## This stuff functions to get the listing of available files
    ## for a given GSE given that there may be many GSEMatrix
    ## files for a given GSE.
    stub = gsub('\\d{1,3}$','nnn',GEO,perl=TRUE)
    gdsurl <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/'
    b = getDirListing(sprintf(gdsurl,stub,GEO))
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
            download.file(sprintf('https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s',
                                  stub,GEO,b[i]),destfile=destfile,mode='wb',
                          method=getOption('download.file.method.GEOquery'))
        }
        ret[[b[i]]] <- parseGSEMatrix(destfile,destdir=destdir,AnnotGPL=AnnotGPL,getGPL=getGPL)$eset
    }
    return(ret)
}

#' Parse a GSE mstrix file
#'
#' Not meant for user calling, but parses
#' a single GSEMatrix file.
#'
#' @importFrom dplyr select filter mutate mutate_all
#' @importFrom tidyr gather spread separate
#' @importFrom readr read_lines
#' @importClassesFrom Biobase ExpressionSet
#' @importFrom magrittr %>%
#'
#' @param fname filename
#' @param AnnotGPL set to TRUE to get the annotation GPL version
#' @param destdir the destdination directory for download
#' @param getGPL whether or not to get the GPL associated
#' @param parseCharacteristics Whether or not to do full "characteristic" parsing
#' @keywords internal
#' 
parseGSEMatrix <- function(fname,AnnotGPL=FALSE,destdir=tempdir(),getGPL=TRUE,parseCharacteristics=TRUE) {
    dat <- read_lines(fname)
    ## get the number of !Series and !Sample lines
    series_header_row_count <- sum(grepl("^!Series_", dat))
    samples_header_row_count <- sum(grepl("^!Sample_", dat))
    series_table_begin_line = grep("^!series_matrix_table_begin", dat)
    if(length(series_table_begin_line) != 1) {
        stop("parsing failed--expected only one '!series_data_table_begin'")
    }
    #con <- fileOpen(fname)
    ## Read the !Series_ and !Sample_ lines
    header <- read.table(fname,sep="\t",header=FALSE,nrows=series_header_row_count)
    tmpdat <- read.table(fname,sep="\t",header=FALSE,nrows=samples_header_row_count,
                         skip=series_header_row_count)
    tmptmp <- t(tmpdat)
    sampledat <- rbind(data.frame(),tmptmp[-1,])
    colnames(sampledat) <- make.unique(sub('!Sample_','',as.character(tmpdat[,1])))
    sampledat[['geo_accession']]=as.character(sampledat[['geo_accession']])
    rownames(sampledat) = sampledat[['geo_accession']]
    ## Lots of GSEs now use "characteristics_ch1" and
    ## "characteristics_ch2" for key-value pairs of
    ## annotation. If that is the case, this simply
    ## cleans those up and transforms the keys to column
    ## names and the values to column values.
    if(length(grep('characteristics_ch',colnames(sampledat)))>0 && parseCharacteristics) {
        pd = sampledat %>%
            dplyr::select(dplyr::contains('characteristics_ch')) %>%
            dplyr::mutate(accession = rownames(.)) %>%
            # these next two lines avoid warnings due
            # to columns having different factor levels
            # (attributes).
            mutate_all(as.character) %>%
            tidyr::gather(characteristics, kvpair, -accession) %>%
            dplyr::filter(grepl(':',kvpair) && !is.na(kvpair))
        # Thx to Mike Smith (@grimbough) for this code
        # sometimes the "characteristics_ch1" fields are empty and contain no 
        # key:value pairs. spread() will fail when called on an
        # empty data_frame.  We catch this case and remove the 
        # "charactics_ch1" column instead
        if(nrow(pd)) {
            pd = dplyr::mutate(pd, characteristics=ifelse(grepl('_ch2',characteristics),'ch2','ch1')) %>%
                tidyr::separate(kvpair, into= c('k','v'), sep=":", fill = 'right', extra = "merge") %>%
                dplyr::mutate(k = paste(k,characteristics,sep=":")) %>%
                dplyr::select(-characteristics) %>%
                dplyr::filter(!is.na(v)) %>%
                tidyr::spread(k, v)
        } else {
            pd = pd %>% 
                dplyr::select(accession)
        }
        ##     dplyr::mutate(characteristics=ifelse(grepl('_ch2',characteristics),'ch2','ch1')) %>%
        ##     dplyr::filter(grepl(':',kvpair)) %>% 
        ##     tidyr::separate(kvpair, into= c('k','v'), sep=":")
        ## if(nrow(pd)>0) {
        ##     pd = pd %>% dplyr::mutate(k = paste(k,characteristics,sep=":")) %>%
        ##         dplyr::select(-characteristics) %>%
        ##         tidyr::spread(k,v)
        sampledat = sampledat %>% dplyr::left_join(pd,by=c('geo_accession'='accession'))
    }
    
    ## used to be able to use colclasses, but some SNP arrays provide only the
    ## genotypes in AA AB BB form, so need to switch it up....
    ##  colClasses <- c('character',rep('numeric',nrow(sampledat)))
    datamat <- read_tsv(fname,quote='"',
                        na=c('NA','null','NULL','Null'), skip = series_table_begin_line,
                        comment = '!series_matrix_table_end')
    tmprownames = datamat[[1]]
                                        # need the as.matrix for single-sample or empty GSE
    datamat <- as.matrix(datamat[!is.na(tmprownames),-1])
    rownames(datamat) <- tmprownames[!is.na(tmprownames)]
    datamat <- as.matrix(datamat)
    rownames(sampledat) <- colnames(datamat)
    GPL=as.character(sampledat[1,grep('platform_id',colnames(sampledat),ignore.case=TRUE)])
    ## if getGPL is FALSE, skip this and featureData is then a data.frame with no columns
    fd = new("AnnotatedDataFrame",data=data.frame(row.names=rownames(datamat)))
    if(getGPL) {
        gpl <- getGEO(GPL,AnnotGPL=AnnotGPL,destdir=destdir)
        vmd <- Columns(gpl)
        dat <- Table(gpl)
        ## GEO uses "TAG" instead of "ID" for SAGE GSE/GPL entries, but it is apparently
        ##     always the first column, so use dat[,1] instead of dat$ID
        ## The next line deals with the empty GSE
        tmpnames=character(0)
        if(ncol(dat)>0) {
          tmpnames=as.character(dat[,1])
        }
        ## Fixed bug caused by an ID being "NA" in GSE15197, for example
        tmpnames[is.na(tmpnames)]="NA"
        rownames(dat) <- make.unique(tmpnames)
        ## Apparently, NCBI GEO uses case-insensitive matching
        ## between platform IDs and series ID Refs ???
        dat <- dat[match(tolower(rownames(datamat)),tolower(rownames(dat))),]
                                        # Fix possibility of duplicate column names in the
                                        # GPL files; this is prevalent in the Annotation GPLs
        rownames(vmd) <- make.unique(colnames(Table(gpl)))
        colnames(dat) <- rownames(vmd)
        fd <- new('AnnotatedDataFrame',data=dat,varMetadata=vmd)
    }
    if(is.null(nrow(datamat))) {
        ## fix empty GSE datamatrix
        ## samplename stuff above does not work with
        ## empty GSEs, so fix here, also
        tmpnames <- names(datamat)
        rownames(sampledat) <- tmpnames
        datamat=matrix(nrow=0,ncol=nrow(sampledat))
        colnames(datamat) <- tmpnames
    } else {
        ## This looks like a dangerous operation but is needed
        ## to deal with the fact that NCBI GEO allows case-insensitive
        ## matching and we need to pick one.
        rownames(datamat) <- rownames(dat)
    }
    eset <- new('ExpressionSet',
                phenoData=as(sampledat,'AnnotatedDataFrame'),
                annotation=GPL,
                featureData=fd,
                exprs=as.matrix(datamat))
    return(list(GPL=as.character(sampledat[1,grep('platform_id',colnames(sampledat),ignore.case=TRUE)]),eset=eset))
}

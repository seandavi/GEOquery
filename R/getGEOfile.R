getGEOfile <- function(GEO,destdir=tempdir(),AnnotGPL=FALSE,
                       amount=c('full','brief','quick','data'))
  {
    amount <- match.arg(amount)
    geotype <- toupper(substr(GEO,1,3))
    mode <- 'wb'
    GEO <- toupper(GEO)
    stub = gsub('\\d{1,3}$','nnn',GEO,perl=TRUE)
    if (geotype == 'GDS') {
      gdsurl <- 'https://ftp.ncbi.nlm.nih.gov/geo/datasets/%s/%s/soft/%s'
      myurl <- sprintf(gdsurl,stub,GEO,paste0(GEO,'.soft.gz'))
      destfile <- file.path(destdir,paste0(GEO,'.soft.gz'))
    }
    if (geotype == 'GSE' & amount=='full') {
      gseurl <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/soft/%s'
      myurl <- sprintf(gseurl,stub,GEO,paste0(GEO,'_family.soft.gz'))
      destfile <- file.path(destdir,paste(GEO,'.soft.gz',sep=""))
    }
    if (geotype == 'GSE' & amount!='full' & amount!='table') {
      gseurl <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
      myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
      destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
      mode <- 'w'
    }
    if (geotype == 'GPL') {
      if (AnnotGPL) {
        gplurl <- 'https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/%s'
        myurl <- sprintf(gplurl,stub,GEO,paste0(GEO,'.annot.gz'))
        destfile <- file.path(destdir,paste(GEO,'.annot.gz',sep=""))
        # check to see if Annotation GPL is present.  If so,
        # use it, else move on to submitter GPL
        res=try({
          if(!file.exists(destfile)) {
            download.file(myurl,destfile,mode=mode,quiet=TRUE,method=getOption('download.file.method.GEOquery'))
            message('File stored at: ')
            message(destfile)
          } else {
            message(sprintf('Using locally cached version of %s found here:\n%s ',GEO,destfile))
          }
        },silent=TRUE)
        if(!inherits(res,'try-error')) {
          return(invisible(destfile))
        } else {
          message('Annotation GPL not available, so will use submitter GPL instead')
        }
      } 
      gseurl <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
      myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
      destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
      mode <- 'w'
      if(!file.exists(destfile)) {
        download.file(myurl,destfile,mode=mode,quiet=TRUE,method=getOption('download.file.method.GEOquery'))
        message('File stored at: ')
        message(destfile)
      } else {
        message(sprintf('Using locally cached version of %s found here:\n%s ',GEO,destfile))
      }
      return(invisible(destfile))
    }
    if (geotype == 'GSM') {
      gseurl <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
      myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
      destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
      mode <- 'w'
    }
    if(!file.exists(destfile)) {
      download.file(myurl,destfile,mode=mode,quiet=TRUE,method=getOption('download.file.method.GEOquery'))
      message('File stored at: ')
      message(destfile)
    } else {
      message(sprintf('Using locally cached version of %s found here:\n%s ',GEO,destfile))
    }
#    if(length(grep('\\.gz',destfile,perl=TRUE))>0) {
#      gunzip(destfile,overwrite=TRUE,remove=TRUE)
#      destfile <- sub('\\.gz$','',destfile)
#    }
    invisible(destfile)
  }

getGEORaw <- function(GEO,destdir=tempdir()) {
  return(getGEOSuppFiles(GEO,baseDir=destdir))
}
                             
gunzip <- function(filename, destname=gsub("[.]gz$", "", filename), overwrite=FALSE, remove=TRUE, BFR.SIZE=1e7) {
  if (filename == destname) 
    stop(sprintf("Argument 'filename' and 'destname' are identical: %s", filename));
  if (!overwrite && file.exists(destname))
    stop(sprintf("File already exists: %s", destname));

  inn <- gzfile(filename, "rb");
  on.exit(if (!is.null(inn)) close(inn));

  out <- file(destname, "wb"); 
  on.exit(close(out), add=TRUE);

  nbytes <- 0;
  repeat { 
    bfr <- readBin(inn, what=raw(0), size=1, n=BFR.SIZE);
    n <- length(bfr);
    if (n == 0)
      break;
    nbytes <- nbytes + n;
    writeBin(bfr, con=out, size=1); 
  };

  if (remove) {
    close(inn);
    inn <- NULL;
    file.remove(filename);
  }
    
  invisible(nbytes);
}

getGEOfile <- function(GEO,destdir=tempdir(),AnnotGPL=FALSE,
                       amount=c('full','brief','quick','data'))
  {
    amount <- match.arg(amount)
    geotype <- toupper(substr(GEO,1,3))
    mode <- 'wb'
    if (geotype == 'GDS') {
      gdsurl <- 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/GDS/'
      myurl <- paste(gdsurl,GEO,'.soft.gz',sep="")
      destfile <- file.path(destdir,paste(GEO,'.soft.gz',sep=""))
    }
    if (geotype == 'GSE' & amount=='full') {
      gseurl <- 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/'
      myurl <- paste(gseurl,GEO,'/',GEO,'_family.soft.gz',sep="")
      destfile <- file.path(destdir,paste(GEO,'.soft.gz',sep=""))
    }
    if (geotype == 'GSE' & amount!='full' & amount!='table') {
      gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
      myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
      destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
      mode <- 'w'
    }
    if (geotype == 'GPL') {
      if (AnnotGPL) {
        gdsurl <- 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/annotation/platforms/'
        myurl <- paste(gdsurl,GEO,'.annot.gz',sep="")
        destfile <- file.path(destdir,paste(GEO,'.annot.gz',sep=""))
      } else {
        gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
        destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
        mode <- 'w'
      }
    }
    if (geotype == 'GSM') {
      gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
      myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
      destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
      mode <- 'w'
    }
    if(!file.exists(destfile)) {
      download.file(myurl,destfile,mode=mode,quiet=TRUE)
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
  geotype <- toupper(substr(GEO,1,3))
  if(geotype=='GSE') {
    GEOurl <- 'ftp://ftp.ncbi.nih.gov/pub/geo/data/geo/raw_data/series/'
    myurl <- paste(GEOurl,GEO,'/',GEO,'_RAW.tar',sep="")
    destfile <- file.path(destdir,paste(GEO,'_RAW.tar',sep=""))
    download.file(myurl,destfile,quiet=TRUE)
    writeLines('File stored at: ')
    writeLines(destfile)
    invisible(destfile)
  } else {
    stop('Fetching raw data supported for GSE only....')
  }
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

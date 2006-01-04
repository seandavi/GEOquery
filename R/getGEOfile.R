getGEOfile <- function(GEO,destdir=tempdir(),
                       amount=c('full','brief','quick','data'))
  {
    amount <- match.arg(amount)
    geotype <- toupper(substr(GEO,1,3))
    mode <- 'wb'
    if (geotype == 'GDS') {
      gdsurl <- 'ftp://ftp.ncbi.nih.gov/pub/geo/data/gds/soft_gz/'
      myurl <- paste(gdsurl,GEO,'.soft.gz',sep="")
      destfile <- file.path(destdir,paste(GEO,'.soft.gz',sep=""))
    }
    if (geotype == 'GSE' & amount=='full') {
      gseurl <- 'ftp://ftp.ncbi.nih.gov/pub/geo/data/geo/by_series/'
      myurl <- paste(gseurl,GEO,'_family.soft.gz',sep="")
      destfile <- file.path(destdir,paste(GEO,'.soft.gz',sep=""))
    }
    if (geotype == 'GSE' & amount!='full' & amount!='table') {
      gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
      myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
      destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
      mode <- 'w'
    }
    if (geotype == 'GPL') {
      gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
      myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
      destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
      mode <- 'w'
    }
    if (geotype == 'GSM') {
      gseurl <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
      myurl <- paste(gseurl,'?targ=self&acc=',GEO,'&form=text&view=',amount,sep='')
      destfile <- file.path(destdir,paste(GEO,'.soft',sep=""))
      mode <- 'w'
    }
    download.file(myurl,destfile,mode=mode)
    writeLines('File stored at: ')
    writeLines(destfile)
    invisible(destfile)
  }

getGEORaw <- function(GEO,destdir=tempdir()) {
  geotype <- toupper(substr(GEO,1,3))
  if(geotype=='GSE') {
    GEOurl <- 'ftp://ftp.ncbi.nih.gov/pub/geo/data/geo/raw_data/series/'
    myurl <- paste(GEOurl,GEO,'/',GEO,'_RAW.tar',sep="")
    destfile <- file.path(destdir,paste(GEO,'_RAW.tar',sep=""))
    download.file(myurl,destfile)
    writeLines('File stored at: ')
    writeLines(destfile)
    invisible(destfile)
  } else {
    stop('Fetching raw data supported for GSE only....')
  }
}

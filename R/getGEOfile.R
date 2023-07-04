#' Download a file from GEO soft file to the local machine
#' 
#' This function simply downloads a SOFT format file associated with the GEO
#' accession number given.
#' 
#' This function downloads GEO SOFT files based on accession number.  It does
#' not do any parsing.  The first two arguments should be fairly
#' self-explanatory, but the last is based on the input to the acc.cgi url at
#' the geo website.  In the default 'full' mode, the entire SOFT format file is
#' downloaded.  Both 'brief' and 'quick' offer shortened versions of the files,
#' good for 'peeking' at the file before a big download on a slow connection.
#' Finally, 'data' downloads only the data table part of the SOFT file and is
#' good for downloading a simple EXCEL-like file for use with other programs (a
#' convenience).
#' 
#' @param GEO Character string, the GEO accession for download (eg., GDS84,
#' GPL96, GSE2553, or GSM10)
#' @param destdir Directory in which to store the resulting downloaded file.
#' Defaults to tempdir()
#' @param AnnotGPL A boolean defaulting to FALSE as to whether or not to use
#' the Annotation GPL information.  These files are nice to use because they
#' contain up-to-date information remapped from Entrez Gene on a regular basis.
#' However, they do not exist for all GPLs; in general, they are only available
#' for GPLs referenced by a GDS
#' @param amount Amount of information to pull from GEO.  Only applies to GSE,
#' GPL, or GSM.  See details...
#' @return Invisibly returns the full path of the downloaded file.
#'
#' @importFrom utils download.file
#'
#' 
#' @author Sean Davis
#' @seealso \code{\link{getGEO}}
#' @references http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
#' @keywords IO
#' @examples
#' 
#'   # myfile <- getGEOfile('GDS10')
#'
#' @export
getGEOfile <- function(GEO, destdir = tempdir(), AnnotGPL = FALSE, amount = c("full",
    "brief", "quick", "data")) {
    amount <- match.arg(amount)
    geotype <- toupper(substr(GEO, 1, 3))
    mode <- "wb"
    GEO <- toupper(GEO)
    stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
    if (geotype == "GDS") {
        gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/datasets/%s/%s/soft/%s"
        myurl <- sprintf(gdsurl, stub, GEO, paste0(GEO, ".soft.gz"))
        destfile <- file.path(destdir, paste0(GEO, ".soft.gz"))
    }
    if (geotype == "GSE" & amount == "full") {
        gseurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/soft/%s"
        myurl <- sprintf(gseurl, stub, GEO, paste0(GEO, "_family.soft.gz"))
        destfile <- file.path(destdir, paste(GEO, ".soft.gz", sep = ""))
    }
    if (geotype == "GSE" & amount != "full" & amount != "table") {
        gseurl <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        myurl <- paste(gseurl, "?targ=self&acc=", GEO, "&form=text&view=", amount,
            sep = "")
        destfile <- file.path(destdir, paste(GEO, ".soft", sep = ""))
        mode <- "w"
    }
    if (geotype == "GPL") {
        if (AnnotGPL) {
            gplurl <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/%s"
            myurl <- sprintf(gplurl, stub, GEO, paste0(GEO, ".annot.gz"))
            destfile <- file.path(destdir, paste(GEO, ".annot.gz", sep = ""))
            # check to see if Annotation GPL is present.  If so, use it, else
            # move on to submitter GPL
            res = try({
                if (!file.exists(destfile)) {
                  downloadFile(myurl, destfile, mode)
                } else {
                  message(sprintf("Using locally cached version of %s found here:\n%s ",
                    GEO, destfile))
                }
            }, silent = TRUE)
            if (!inherits(res, "try-error")) {
                return(invisible(destfile))
            } else {
                message("Annotation GPL not available, so will use submitter GPL instead")
            }
        }
        gseurl <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        myurl <- paste(gseurl, "?targ=self&acc=", GEO, "&form=text&view=", amount,
            sep = "")
        destfile <- file.path(destdir, paste(GEO, ".soft.gz", sep = ""))
        mode <- "w"
        if (!file.exists(destfile)) {
            downloadFile(myurl, destfile, mode = mode, quiet = TRUE)
        } else {
            message(sprintf("Using locally cached version of %s found here:\n%s ",
                GEO, destfile))
        }
        return(invisible(destfile))
    }
    if (geotype == "GSM") {
        gseurl <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        myurl <- paste(gseurl, "?targ=self&acc=", GEO, "&form=text&view=", amount,
            sep = "")
        destfile <- file.path(destdir, paste(GEO, ".soft", sep = ""))
        mode <- "w"
    }
    if (!file.exists(destfile)) {
        downloadFile(myurl, destfile, mode)
    } else {
        message(sprintf("Using locally cached version of %s found here:\n%s ", GEO,
            destfile))
    }
    # if(length(grep('\\.gz',destfile,perl=TRUE))>0) {
    # gunzip(destfile,overwrite=TRUE,remove=TRUE) destfile <-
    # sub('\\.gz$','',destfile) }
    invisible(destfile)
}

getGEORaw <- function(GEO, destdir = tempdir()) {
    return(getGEOSuppFiles(GEO, baseDir = destdir))
}



#' Gunzip a file
#' 
#' gunzip a file
#' 
#' This function was stripped out of R.utils due to breaking some stuff on the
#' bioconductor build machine.
#' 
#' @param filename The filename to be unzipped
#' @param destname The destination file
#' @param overwrite Boolean indicating whether or not to overwrite a destfile
#' of the same name
#' @param remove Boolean indicating whether or not to remove the original file
#' after completion
#' @param BFR.SIZE The size of the read buffer....
#' @return Invisibly, the number of bytes read.
#' @author Original author: Henrik Bengtsson
#' @seealso \code{\link{gzfile}}
#' @keywords IO
#'
#' @examples
#' 
#' # gunzip('file.gz',remove=FALSE)
#' 
#' 
#' @export
gunzip <- function(filename, destname = gsub("[.]gz$", "", filename), overwrite = FALSE,
    remove = TRUE, BFR.SIZE = 1e+07) {
    if (filename == destname)
        stop(sprintf("Argument 'filename' and 'destname' are identical: %s", filename))
    if (!overwrite && file.exists(destname))
        stop(sprintf("File already exists: %s", destname))

    inn <- gzfile(filename, "rb")
    on.exit(if (!is.null(inn)) close(inn))

    out <- file(destname, "wb")
    on.exit(close(out), add = TRUE)

    nbytes <- 0
    repeat {
        bfr <- readBin(inn, what = raw(0), size = 1, n = BFR.SIZE)
        n <- length(bfr)
        if (n == 0)
            break
        nbytes <- nbytes + n
        writeBin(bfr, con = out, size = 1)
    }

    if (remove) {
        close(inn)
        inn <- NULL
        file.remove(filename)
    }

    invisible(nbytes)
}

# internal use only
downloadFile <- function(url, destfile, mode, quiet = TRUE) {
    h <- curl::new_handle()
    curl::handle_setheaders(h, `accept-encoding` = "gzip")
    timeout_seconds <- max(getOption("timeout"), 120)
    curl::handle_setopt(h, timeout_ms = timeout_seconds * 1000)
    result = tryCatch({
        curl::curl_download(url, destfile, mode = mode, quiet = quiet, handle = h)
        return(TRUE)
    }, error = function(e) {
        message(e)
        return(FALSE)
    })
    message("File stored at:")
    message(destfile)

    ## if the download failed, remove the corrupted file and report the error
    if (!result) {
        if (file.exists(destfile)) {
            file.remove(destfile)
        }
        stop(sprintf("Failed to download %s!", destfile))
    }
    return(0)
}

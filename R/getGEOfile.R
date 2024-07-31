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
#' # myfile <- getGEOfile('GDS10')
#'
#' @export
getGEOfile <- function(GEO, destdir = tempdir(), AnnotGPL = FALSE, amount = c(
                         "full",
                         "brief", "quick", "data"
                       )) {
  amount <- match.arg(amount)
  geotype <- toupper(substr(GEO, 1, 3))
  mode <- "wb"
  GEO <- toupper(GEO)
  stub <- gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
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
      sep = ""
    )
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
      res <- try(
        {
          if (!file.exists(destfile)) {
            downloadFile(myurl, destfile, mode)
          } else {
            message(sprintf(
              "Using locally cached version of %s found here:\n%s ",
              GEO, destfile
            ))
          }
        },
        silent = TRUE
      )
      if (!inherits(res, "try-error")) {
        return(invisible(destfile))
      } else {
        message("Annotation GPL not available, so will use submitter GPL instead")
      }
    }
    gseurl <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    myurl <- paste(gseurl, "?targ=self&acc=", GEO, "&form=text&view=", amount,
      sep = ""
    )
    destfile <- file.path(destdir, paste(GEO, ".soft.gz", sep = ""))
    mode <- "w"
    if (!file.exists(destfile)) {
      downloadFile(myurl, destfile, mode = mode, quiet = TRUE)
    } else {
      message(sprintf(
        "Using locally cached version of %s found here:\n%s ",
        GEO, destfile
      ))
    }
    return(invisible(destfile))
  }
  if (geotype == "GSM") {
    gseurl <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    myurl <- paste(gseurl, "?targ=self&acc=", GEO, "&form=text&view=", amount,
      sep = ""
    )
    destfile <- file.path(destdir, paste(GEO, ".soft", sep = ""))
    mode <- "w"
  }
  if (!file.exists(destfile)) {
    downloadFile(myurl, destfile, mode)
  } else {
    message(sprintf(
      "Using locally cached version of %s found here:\n%s ", GEO,
      destfile
    ))
  }
  invisible(destfile)
}

getGEORaw <- function(GEO, destdir = tempdir()) {
  return(getGEOSuppFiles(GEO, baseDir = destdir))
}


# internal use only
downloadFile <- function(url, destfile, mode, quiet = TRUE) {
  timeout_seconds <- max(getOption("timeout"), 120)

  result <- tryCatch(
    {
      req <- httr2::request(url) %>%
        httr2::req_headers(`accept-encoding` = "gzip") %>%
        httr2::req_timeout(timeout_seconds)

      resp <- httr2::req_perform(req)

      if (httr2::resp_status(resp) == 200) {
        writeBin(httr2::resp_body_raw(resp), destfile)
        return(TRUE)
      } else {
        stop("Failed to download file: ", httr2::resp_status(resp))
      }
    },
    error = function(e) {
      message(e)
      return(FALSE)
    }
  )

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

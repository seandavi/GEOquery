# Join a base URL and a path (or vector of paths) with a single slash,
# preserving the scheme's '//'. Unlike file.path(), this does not collapse
# 'https://' into 'https:/' and does not produce a double slash when the base
# already ends in '/' (#131). Vectorized over `path`.
url_join <- function(base, path) {
    paste0(sub("/+$", "", base), "/", sub("^/+", "", path))
}

#' get a directory listing from NCBI GEO
#'
#' This one makes some assumptions about the structure of the HTML response
#' returned.
#'
#' @param url A URL, assumed to return an NCBI-formatted
#'   index page
#' 
#' @importFrom xml2 read_html xml_text xml_find_all 
getDirListing <- function(url) {
    # Fetch the index page through the shared httr2 request (.geo_request) so it
    # uses the same timeout/retry handling as downloadFile and can be mocked in
    # tests (#173). Returns a character vector of GEO filenames.
    resp <- httr2::req_perform(.geo_request(url))
    a <- xml2::read_html(httr2::resp_body_string(resp))
    fnames <- grep("^G", xml2::xml_text(xml2::xml_find_all(a, "//a/@href")), value = TRUE)
    return(fnames)
}

#' Get GEO supplemental file URL for a given GEO accession
#'
#' @param GEO
#'
#' @examples
#' # an example of a GEO supplemental file URL
#' # with a set of single-cell RNA-seq data
#' url = getGEOSuppFileURL("GSE161228")
#' url
#' 
#' \dontrun{
#'   browseURL(url)
#' }
#' 
#' @export
getGEOSuppFileURL <- function(GEO) {
    stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
    geotype <- toupper(substr(GEO, 1, 3))
    if (geotype == "GSM") {
        url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/", stub,
            GEO)
    }
    if (geotype == "GSE") {
        url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", stub,
            GEO)
    }
    if (geotype == "GPL") {
        url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platform/%s/%s/suppl/", stub,
            GEO)
    }
    return(url)
}



#' Get Supplemental Files from GEO
#' 
#' NCBI GEO allows supplemental files to be attached to GEO Series (GSE), GEO
#' platforms (GPL), and GEO samples (GSM).  This function 'knows' how to get
#' these files based on the GEO accession.  No parsing of the downloaded files
#' is attempted, since the file format is not generally knowable by the
#' computer.
#' 
#' Again, just a note that the files are simply downloaded.
#' 
#' @param GEO A GEO accession number such as GPL1073 or GSM1137
#' @param makeDirectory Should a 'subdirectory' for the downloaded files be
#' created?  Default is TRUE.  If FALSE, the files will be downloaded directly
#' into the baseDir.
#' @param baseDir The base directory for the downloads.  Default is the current
#' working directory.
#' @param filter_regex A character(1) regular expression that will be
#'     used to filter the filenames from GEO to limit those files that
#'     will be downloaded. This is useful to limit to, for example,
#'     bed files only.
#' @param fetch_files logical(1). If TRUE, then actually download the
#'     files. If FALSE, just return the filenames that would have been
#'     downloaded. Useful for testing and getting a list of files
#'     without actual download.
#' @param quiet logical(1). If TRUE, suppress informational messages such as
#'     "No supplemental files found" and "Using locally cached version".
#'     Defaults to the `GEOquery.quiet` option, or FALSE.
#' @return If fetch_files=TRUE, a data frame is returned invisibly with rownames representing the
#' full path of the resulting downloaded files and the records in the
#' data.frame the output of file.info for each downloaded file.
#' If fetch_files=FALSE, a data.frame of URLs and filenames is returned.
#' @author Sean Davis <sdavis2@@mail.nih.gov>
#' @keywords IO database
#' @examples
#' \dontrun{
#' 
#' a <- getGEOSuppFiles('GSM1137', fetch_files = FALSE)
#' a
#' 
#' # with a set of single-cell RNA-seq data
#' a <- getGEOSuppFiles('GSE161228', fetch_files = FALSE)
#' a
#' 
#' }
#' @export
getGEOSuppFiles <- function(
    GEO,
    makeDirectory = TRUE,
    baseDir = getwd(),
    fetch_files = TRUE,
    filter_regex = NULL,
    quiet = getOption("GEOquery.quiet", FALSE)
) {
    geotype <- toupper(substr(GEO, 1, 3))
    storedir <- baseDir
    fileinfo <- list()
    url <- getGEOSuppFileURL(GEO)
    fnames <- try(getDirListing(url), silent = TRUE)
    if (inherits(fnames, "try-error")) {
        if (!quiet) {
            message("No supplemental files found.")
            message("Check URL manually if in doubt")
            message(url)
        }
        return(NULL)
    }
    if (makeDirectory) {
        suppressWarnings(dir.create(storedir <- file.path(baseDir, GEO)))
    }
    if (!is.null(filter_regex)) {
        fnames = fnames[grepl(filter_regex, fnames)]
    }
    if (fetch_files) {
        for (i in fnames) {
            destfile <- file.path(storedir, i)


 
                if (!file.exists(destfile)) {
                    httr2::request(base_url=url) |>
                      httr2::req_url_path_append(i) |>
                      httr2::req_url_query(tool="geoquery") |>
                      httr2::req_perform(path=destfile)
                    res <- 0
                } else {
                  if (!quiet) {
                    message(sprintf("Using locally cached version of supplementary file(s) %s found here:\n%s ",
                      GEO, destfile))
                  }
                    res <- 0
                }
            fileinfo[[destfile]] <- file.info(destfile)
        }
        ret <- do.call(rbind, fileinfo)
        ret$fname <- fnames
        ret$destdir <- storedir
        ret$filepath <- file.path(storedir, fnames)
        ret$GEO <- GEO
        return(ret)
    } else {
        return(data.frame(fname = fnames, url = url_join(url, fnames)))
    }
}

#' GSE Supplemental file listing
#' 
#' The GEO Series records often have one or more supplemental files.
#' In most cases, those files are archived as '.tar' files, the contents
#' of which are only available in a file listing file not present on the
#' website for download.
#' 
#' This function reads that file listing file and returns the results
#' as a data.frame. 
#' 
#' @returns A data.frame with 5 columns. See example. 
#' 
#' @param GSE character(1) the GSE accession
#' 
#' @examples
#' \dontrun{
#' getGEOSeriesFileListing('GSE288770')
#' 
#' }
#' @export
getGEOSeriesFileListing <- function(GSE) {
  url = getGEOSuppFileURL(GSE)
  ret <- readr::read_tsv(url_join(url, 'filelist.txt'))
  ret |> 
    dplyr::rename_with(tolower) |>
    dplyr::rename('archive_or_file'='#archive/file')
}

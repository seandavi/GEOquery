#' get all download links from a GEO accession
#'
#' This function gets all download links from a GEO accession number.
#'
#' @param gse GEO accession number
#'
#' @return A character vector with all download links
#'
#' @keywords internal
getGSEDownloadPageURLs <- function(gse) {
  url <- "https://ncbi.nlm.nih.gov/geo/download/"
  links <- httr2::request(url) |>
    httr2::req_timeout(15) |>
    httr2::req_retry(3) |>
    httr2::req_url_query(acc = gse) |>
    httr2::req_perform() |>
    httr2::resp_body_string() |>
    rvest::read_html() |>
    rvest::html_nodes("a") |>
    rvest::html_attr("href") |>
    stringr::str_replace("^/geo/", "https://www.ncbi.nlm.nih.gov/geo/") |>
    stringr::str_replace("^ftp://", "https://")
  class(links) <- c("geoDownloadLinks", class(links))
  return(links)
}


#' Get the link to the raw counts file from GEO
#'
#' This function extracts the link to the raw counts file from a
#' geoDownloadLinks object.
#'
#' @param links A geoDownloadLinks object
#'
#' @return A character vector with the link to the raw counts file
#'
#' @keywords internal
getRNAQuantRawCountsURL <- function(links) {
  if (!inherits(links, "geoDownloadLinks")) {
    stop("Input must be a geoDownloadLinks object")
  }
  link <- stringr::str_subset(links, "raw_counts")
  return(link)
}

#' Get the RNA-seq quantification annotation link
#'
#' This function extracts the link to the RNA-seq quantification annotation
#' file from a geoDownloadLinks object.
#'
#' @param links A geoDownloadLinks object
#'
#' @return A character vector with the link to the annotation file
#'
#' @keywords internal
getRNAQuantAnnotationURL <- function(links) {
  if (!inherits(links, "geoDownloadLinks")) {
    stop("Input must be a geoDownloadLinks object")
  }
  link <- stringr::str_subset(links, "annot.tsv.gz")
  return(link)
}


#' Read RNA-seq quantification annotation from GEO
#'
#' This function reads the annotation file from a GEO link. The annotation
#' file is expected to be a tab-separated file with the first column
#' containing the gene IDs and the remaining columns containing the
#' annotation information.
#'
#' @param link A link to the annotation file
#'
#' @return A data frame of annotation information with gene IDs as row names
#'
#' @keywords internal
readRNAQuantAnnotation <- function(link) {
  annotation <- as.data.frame(readr::read_tsv(link, show_col_types = FALSE))
  rownames(annotation) <- as.character(annotation$GeneID)
  return(annotation)
}


#' Read raw counts from GEO
#'
#' This function reads the raw counts from a GEO link. The raw counts
#' are expected to be in a tab-separated file with the first column
#' containing the gene IDs and the remaining columns containing the
#' raw counts.
#'
#' This function reads the raw counts and returns a matrix with the
#' gene IDs as the row names, ready for use in creating a
#' SummarizedExperiment.
#'
#' @param link A link to the raw counts file
#'
#' @return A matrix of raw counts with gene IDs as row names
#'
#' @keywords internal
readRNAQuantRawCounts <- function(link) {
  quants <- readr::read_tsv(link, show_col_types = FALSE)
  gene_ids <- as.character(quants$GeneID)
  quants <- as.matrix(quants[, -1])
  rownames(quants) <- gene_ids
  return(quants)
}


#' Get RNA-seq quantification and annotation from GEO
#'
#' This function downloads the raw counts and annotation files from GEO
#' for a given GEO accession number.
#'
#' @param gse GEO accession number
#'
#' @return A list with two elements: quants (a matrix of raw counts) and
#' annotation (a data frame of annotation information).
#'
#' @keywords internal
getRNASeqQuantResults <- function(gse) {
  links <- getGSEDownloadPageURLs(gse)
  raw_counts_link <- getRNAQuantRawCountsURL(links)
  if (length(raw_counts_link) == 0) {
    stop(
      "No raw counts file found.\n",
      "Navigate to: \n  https://ncbi.nlm.nih.gov/geo/download/?acc=",
      gse,
      "\nand check if the 'RNA-Seq raw counts' link is available."
    )
  }
  annotation_link <- getRNAQuantAnnotationURL(links)
  quants <- readRNAQuantRawCounts(raw_counts_link)
  annotation <- readRNAQuantAnnotation(annotation_link)
  return(list(quants = quants, annotation = annotation))
}

#' Browse GEO search website for RNA-seq datasets
#'
#' This function opens a browser window to the NCBI GEO website
#' with a search for RNA-seq datasets. It is included as a convenience
#' function to remind users of how to search for RNA-seq datasets using
#' the NCBI GEO website and an "rnaseq counts" filter.
#'
#' @examples
#' \dontrun{
#' browseWebsiteRNASeqSearch()
#' }
#' @export
browseWebsiteRNASeqSearch <- function() {
  browseURL(
    "https://ncbi.nlm.nih.gov/gds?term=%22rnaseq%20counts%22%5BFilter%5D"
  )
}

#' Get GEO RNA-seq quantifications as a SummarizedExperiment object
#'
#' For human and mouse GEO datasets, NCBI GEO attempts to process
#' the raw data and provide quantifications in the form of raw counts
#' and an annotation file. This function downloads the raw counts and
#' annotation files from GEO and merges that with the metadata from the GEO
#' object to create a SummarizedExperiment.
#'
#' @details
#' A major barrier to fully exploiting and reanalyzing the massive volumes
#' of public RNA-seq data archived by SRA is the cost and effort required to
#' consistently process raw RNA-seq reads into concise formats that summarize
#' the expression results. To help address this need, the NCBI SRA and GEO
#' teams have built a pipeline that precomputes RNA-seq gene expression counts
#' and delivers them as count matrices that may be incorporated into commonly
#' used differential expression analysis and visualization software.
#'
#' The pipeline processes RNA-seq data from SRA using the HISAT2 aligner and
#' and then generates gene expression counts using the featureCounts program.
#'
#' See the
#' [GEO documentation](https://ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html)
#' for more details.
#'
#' @import SummarizedExperiment
#'
#'
#' @param accession GEO accession number
#'
#' @return A SummarizedExperiment object with the raw counts as the counts
#' assay, the annotation as the rowData, and the metadata from GEO as
#' the colData.
#'
#'
#' @examples
#' se <- getRNASeqData("GSE164073")
#' se
#'
#' @export
getRNASeqData <- function(accession) {
  quantifications <- getRNASeqQuantResults(accession)
  se <- as(
    getGEO(accession)[[1]],
    "SummarizedExperiment"
  )
  new_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = quantifications$quants),
    rowData = quantifications$annotation,
    colData = SummarizedExperiment::colData(se),
    metadata = S4Vectors::metadata(se)
  )
  new_se
}

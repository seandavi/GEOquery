# Single-cell support (#158, ADR-0004).
#
# Single-cell data on GEO lives in supplementary files, not the Series Matrix,
# and the file naming is highly varied. geoSingleCellManifest() inventories a
# GSE's supplementary files and classifies them by single-cell format so a user
# can see what a study contains -- and how 10x triplets group by sample --
# before downloading what may be many gigabytes. The classification logic is
# factored into pure helpers so it can be tested without network access.

# Classify a single supplementary filename into a (format, role) pair.
.classify_sc_file <- function(fname) {
    low <- tolower(fname)
    if (grepl("matrix\\.mtx(\\.gz)?$", low)) return(c(format = "10x_mtx", role = "matrix"))
    if (grepl("barcodes\\.tsv(\\.gz)?$", low)) return(c(format = "10x_mtx", role = "barcodes"))
    if (grepl("(features|genes)\\.tsv(\\.gz)?$", low)) return(c(format = "10x_mtx", role = "features"))
    if (grepl("\\.h5ad$", low)) return(c(format = "h5ad", role = "anndata"))
    if (grepl("\\.h5$", low)) return(c(format = "10x_h5", role = "matrix"))
    if (grepl("\\.loom$", low)) return(c(format = "loom", role = "matrix"))
    if (grepl("\\.(rds|rdata)$", low)) return(c(format = "rds", role = "object"))
    if (grepl("\\.tar(\\.gz)?$", low)) return(c(format = "tar", role = "archive"))
    c(format = "other", role = NA_character_)
}

# Extract the GSM sample id from a filename, if present.
.sc_sample_id <- function(fname) {
    m <- regmatches(fname, regexpr("GSM[0-9]+", fname))
    if (length(m) == 0) NA_character_ else m[1]
}

# Build the manifest data.frame from a vector of filenames (and optional URLs).
.classify_sc_files <- function(fnames, urls = NA_character_) {
    cls <- t(vapply(fnames, .classify_sc_file, character(2)))
    data.frame(
        fname = as.character(fnames),
        sample = vapply(fnames, .sc_sample_id, character(1)),
        format = cls[, "format"],
        role = cls[, "role"],
        url = urls,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

#' Inventory the single-cell supplementary files of a GEO Series
#'
#' Lists the supplementary files attached to a GSE and classifies each by
#' single-cell format (10x Matrix Market triplet, 10x HDF5, AnnData h5ad, loom,
#' Seurat rds, tar archive, or other), extracting the GSM sample id where
#' present. This lets you see what a single-cell study contains -- and how 10x
#' triplets group by sample -- before downloading potentially many gigabytes.
#'
#' No files are downloaded. The result feeds the planned single-cell readers
#' (see ADR-0004); reading itself uses Bioconductor importers (TENxIO, anndataR)
#' that are optional dependencies.
#'
#' @param GEO A GEO Series accession, e.g. "GSE161228".
#' @return A data.frame with columns \code{fname}, \code{sample} (GSM id or NA),
#'   \code{format}, \code{role}, and \code{url}. Zero rows if the GSE has no
#'   supplementary files.
#' @seealso \code{\link{getGEOSuppFiles}}
#' @examples
#' \dontrun{
#'   m <- geoSingleCellManifest("GSE161228")
#'   m
#' }
#' @export
geoSingleCellManifest <- function(GEO) {
    files <- getGEOSuppFiles(GEO, fetch_files = FALSE, quiet = TRUE)
    if (is.null(files) || nrow(files) == 0) {
        return(.classify_sc_files(character(0), character(0)))
    }
    .classify_sc_files(files$fname, files$url)
}

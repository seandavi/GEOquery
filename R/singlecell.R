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
    if (length(fnames) == 0) {
        return(data.frame(
            fname = character(0), sample = character(0), format = character(0),
            role = character(0), url = character(0), stringsAsFactors = FALSE
        ))
    }
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

# Roles that make a complete 10x Matrix Market unit.
.mtx_roles <- c("matrix", "barcodes", "features")

# Group a manifest into loadable units (one per sample + format) and assess
# completeness. Pure logic; used by geoSingleCellUnits().
.sc_units <- function(manifest) {
    cols <- c("sample", "format", "n_files", "status", "loadable")
    if (nrow(manifest) == 0) {
        out <- data.frame(matrix(nrow = 0, ncol = length(cols)))
        colnames(out) <- cols
        return(out)
    }
    key <- paste(manifest$sample, manifest$format, sep = "|")
    parts <- lapply(split(manifest, key), function(g) {
        fmt <- g$format[1]
        if (fmt == "10x_mtx") {
            missing <- setdiff(.mtx_roles, unique(g$role))
            status <- if (length(missing) == 0) {
                "complete"
            } else {
                sprintf("incomplete (missing %s)", paste(missing, collapse = ", "))
            }
        } else if (fmt %in% c("h5ad", "10x_h5", "loom", "rds")) {
            status <- "complete"
        } else {
            status <- "unsupported"
        }
        data.frame(
            sample = g$sample[1], format = fmt, n_files = nrow(g),
            status = status, loadable = identical(status, "complete"),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out[order(out$sample, out$format), ]
}

#' Group a single-cell manifest into loadable units
#'
#' Collapses a \code{\link{geoSingleCellManifest}} into one row per loadable
#' unit (a sample + format combination) and reports completeness. A 10x Matrix
#' Market unit is "complete" only when its matrix, barcodes, and features files
#' are all present; single-file formats (h5ad, 10x h5, loom, rds) are always
#' complete. The \code{loadable} column flags units a reader can consume.
#'
#' @param manifest A data.frame returned by \code{geoSingleCellManifest()}.
#' @return A data.frame with columns \code{sample}, \code{format},
#'   \code{n_files}, \code{status}, and \code{loadable}.
#' @seealso \code{\link{geoSingleCellManifest}}
#' @examples
#' \dontrun{
#'   m <- geoSingleCellManifest("GSE161228")
#'   geoSingleCellUnits(m)
#' }
#' @export
geoSingleCellUnits <- function(manifest) {
    .sc_units(manifest)
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


# ---- Readers (#158, SC3; ADR-0004) ----------------------------------------
#
# Reading uses Bioconductor importers, kept as optional (Suggests) dependencies
# behind requireNamespace() guards:
#   10x Matrix Market / 10x HDF5 -> TENxIO
#   AnnData (.h5ad)              -> anndataR
# NOT covered (read with their native packages): loom, Seurat .rds, and
# idiosyncratic layouts (e.g. a single combined matrix for many samples).

.require_pkg <- function(pkg, what) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf(
            "Reading %s requires the '%s' package. Install it with BiocManager::install('%s').",
            what, pkg, pkg
        ), call. = FALSE)
    }
}

# Arrange a 10x triplet into a directory with canonical filenames that
# TENxIO::TENxFileList() expects. `x` may already be such a directory.
.arrange_10x <- function(x) {
    if (length(x) == 1 && dir.exists(x)) {
        return(x)
    }
    d <- tempfile("tenx_")
    dir.create(d)
    for (f in x) {
        role <- .classify_sc_file(f)[["role"]]
        gz <- if (grepl("\\.gz$", f)) ".gz" else ""
        canon <- switch(role,
            matrix = paste0("matrix.mtx", gz),
            barcodes = paste0("barcodes.tsv", gz),
            features = paste0("features.tsv", gz),
            basename(f)
        )
        file.copy(f, file.path(d, canon))
    }
    d
}

.read_sc_10x_mtx <- function(x) {
    .require_pkg("TENxIO", "10x Matrix Market data")
    BiocIO::import(TENxIO::TENxFileList(.arrange_10x(x)))
}
.read_sc_10x_h5 <- function(path) {
    .require_pkg("TENxIO", "10x HDF5 data")
    BiocIO::import(TENxIO::TENxH5(path))
}
.read_sc_h5ad <- function(path) {
    .require_pkg("anndataR", "AnnData (h5ad) data")
    anndataR::read_h5ad(path, as = "SingleCellExperiment")
}

# Read a `.rds` supplementary file. Its contents are opaque from the filename,
# so detect by class: a saved SingleCellExperiment is returned as-is; a saved
# Seurat object is coerced to SingleCellExperiment (the package lingua franca)
# via Seurat; anything else is an error (#197).
.read_sc_rds <- function(path) {
    obj <- readRDS(path)
    if (methods::is(obj, "SingleCellExperiment")) {
        return(obj)
    }
    if (inherits(obj, "Seurat")) {
        .require_pkg("Seurat", "Seurat .rds files")
        return(Seurat::as.SingleCellExperiment(obj))
    }
    stop(sprintf(
        paste0("Unsupported .rds contents (class: %s). readGEOSingleCell() ",
            "reads .rds files containing a Seurat or SingleCellExperiment object."),
        paste(class(obj), collapse = "/")
    ), call. = FALSE)
}

# Coerce a SingleCellExperiment to the requested output class (#196). Seurat
# output is produced by coercion and requires the Seurat package.
.as_sc_output <- function(sce, as = c("SingleCellExperiment", "Seurat")) {
    as <- match.arg(as)
    if (as == "Seurat") {
        .require_pkg("Seurat", "Seurat output")
        return(Seurat::as.Seurat(sce))
    }
    sce
}

#' Read a single-cell file (or 10x triplet) into a SingleCellExperiment
#'
#' Low-level reader: given already-downloaded local file(s), dispatch on format
#' to the appropriate Bioconductor importer and return a
#' \code{SingleCellExperiment} (or a \code{Seurat} object with
#' \code{as = "Seurat"}). Use this for full control; see
#' \code{\link{getGEOSingleCell}} for the high-level convenience wrapper.
#'
#' Supported formats: \code{"10x_mtx"} (a directory, or the matrix/barcodes/
#' features files, read via TENxIO), \code{"10x_h5"} (CellRanger HDF5, TENxIO),
#' \code{"h5ad"} (AnnData, anndataR), and \code{"rds"} (a saved
#' \code{SingleCellExperiment} or \code{Seurat} object; detected by class).
#' loom is not supported -- read it with \code{LoomExperiment} directly.
#'
#' @param x A path to a single file (\code{.h5}/\code{.h5ad}/\code{.rds}), a
#'   directory containing a 10x triplet, or a character vector of the triplet
#'   files.
#' @param format One of "10x_mtx", "10x_h5", "h5ad", "rds". If NULL (default),
#'   guessed from \code{x}.
#' @param as Output class, one of "SingleCellExperiment" (default) or "Seurat"
#'   (coerced via the Seurat package, an optional dependency).
#' @return A \code{SingleCellExperiment}, or a \code{Seurat} object if
#'   \code{as = "Seurat"}.
#' @seealso \code{\link{getGEOSingleCell}}, \code{\link{geoSingleCellManifest}}
#' @export
readGEOSingleCell <- function(x, format = NULL,
    as = c("SingleCellExperiment", "Seurat")) {
    as <- match.arg(as)
    if (is.null(format)) {
        format <- .classify_sc_file(x[1])[["format"]]
    }
    sce <- switch(format,
        "10x_mtx" = .read_sc_10x_mtx(x),
        "10x_h5" = .read_sc_10x_h5(x),
        "h5ad" = .read_sc_h5ad(x),
        "rds" = .read_sc_rds(x),
        stop(sprintf(
            paste0("Single-cell format '%s' is not supported by readGEOSingleCell(). ",
                "Supported: 10x_mtx, 10x_h5, h5ad, rds. loom is not handled; ",
                "read it with LoomExperiment directly."),
            format
        ), call. = FALSE)
    )
    .as_sc_output(sce, as)
}

# When a sample offers more than one loadable format, keep just one, by
# preference (richer/standard first).
.format_priority <- c("h5ad", "10x_h5", "10x_mtx")
.prefer_one_format <- function(load) {
    if (nrow(load) == 0) {
        return(load)
    }
    keep <- lapply(split(load, load$sample), function(g) {
        if (nrow(g) == 1) {
            return(g)
        }
        ord <- order(match(g$format, .format_priority))
        g[ord[1], , drop = FALSE]
    })
    out <- do.call(rbind, keep)
    rownames(out) <- NULL
    out[order(out$sample), , drop = FALSE]
}

# Pure unit-selection logic: split units into those to load and those skipped,
# honoring optional `samples` / `format` filters and one-format-per-sample.
.select_sc_units <- function(units, samples = NULL, format = NULL) {
    load <- units[units$loadable %in% TRUE, , drop = FALSE]
    if (!is.null(samples)) {
        load <- load[load$sample %in% samples, , drop = FALSE]
    }
    if (!is.null(format)) {
        load <- load[load$format %in% format, , drop = FALSE]
    }
    load <- .prefer_one_format(load)
    loaded_key <- paste(load$sample, load$format)
    skip <- units[!(paste(units$sample, units$format) %in% loaded_key), , drop = FALSE]
    list(load = load, skip = skip)
}

#' Download and read the single-cell data of a GEO Series
#'
#' High-level, best-effort convenience wrapper: inventories the GSE
#' (\code{\link{geoSingleCellManifest}}), groups files into loadable units
#' (\code{\link{geoSingleCellUnits}}), downloads each loadable unit, reads it
#' with \code{\link{readGEOSingleCell}}, and returns the results. It reports
#' which units it loads and which it skips.
#'
#' This handles common, well-structured layouts (clean per-sample 10x, h5ad, or
#' a saved object in \code{.rds}). It does NOT handle every GSE: loom files,
#' files packaged inside a \code{_RAW.tar} archive, and idiosyncratic layouts
#' (e.g. a single combined matrix for many samples) are out of scope -- use the
#' manifest plus \code{readGEOSingleCell()} directly for those.
#'
#' @param GEO A GEO Series accession, e.g. "GSE161228".
#' @param samples Optional character vector of GSM ids to restrict to.
#' @param format Optional format(s) to restrict to ("10x_mtx", "10x_h5",
#'   "h5ad", "rds").
#' @param combine Logical; if TRUE attempt to \code{cbind} the per-sample
#'   objects into one (requires matching features). Default FALSE returns a list.
#' @param as Output class, one of "SingleCellExperiment" (default) or "Seurat"
#'   (coerced via the Seurat package, an optional dependency).
#' @param destdir Download destination directory.
#' @return A named list of \code{SingleCellExperiment} (one per sample), or a
#'   single combined object if \code{combine = TRUE}; \code{Seurat} objects if
#'   \code{as = "Seurat"}.
#' @seealso \code{\link{geoSingleCellManifest}}, \code{\link{readGEOSingleCell}}
#' @export
getGEOSingleCell <- function(GEO, samples = NULL, format = NULL, combine = FALSE,
    as = c("SingleCellExperiment", "Seurat"), destdir = tempdir()) {
    as <- match.arg(as)
    manifest <- geoSingleCellManifest(GEO)
    units <- geoSingleCellUnits(manifest)
    sel <- .select_sc_units(units, samples, format)
    if (nrow(sel$load) == 0) {
        stop(sprintf(
            "No loadable single-cell units found for %s. Inspect geoSingleCellManifest('%s').",
            GEO, GEO
        ), call. = FALSE)
    }
    if (nrow(sel$skip) > 0) {
        message(sprintf(
            "Skipping %d unit(s): %s", nrow(sel$skip),
            paste(sprintf("%s [%s]", sel$skip$sample, sel$skip$status), collapse = "; ")
        ))
    }
    results <- list()
    for (i in seq_len(nrow(sel$load))) {
        u <- sel$load[i, ]
        message(sprintf("Loading %s (%s)...", u$sample, u$format))
        unit_files <- manifest[manifest$sample %in% u$sample & manifest$format == u$format, ]
        dl <- getGEOSuppFiles(GEO, fetch_files = TRUE, baseDir = destdir,
            filter_regex = u$sample, quiet = TRUE)
        local <- dl$filepath[basename(dl$filepath) %in% unit_files$fname]
        # read as SingleCellExperiment (the lingua franca); coerce on output
        results[[u$sample]] <- readGEOSingleCell(local, format = u$format)
    }
    if (combine && length(results) > 1) {
        return(.as_sc_output(do.call(SummarizedExperiment::cbind, results), as))
    }
    lapply(results, .as_sc_output, as = as)
}

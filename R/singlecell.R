# Single-cell support (#158, ADR-0004).
#
# Single-cell data on GEO lives in supplementary files, not the Series Matrix,
# and the file naming is highly varied. geoSingleCellManifest() inventories a
# GSE's supplementary files and classifies them by single-cell format so a user
# can see what a study contains -- and how 10x triplets group by sample --
# before downloading what may be many gigabytes. The classification logic is
# factored into pure helpers so it can be tested without network access.
#
# A GSE's series-level suppl directory often does NOT hold the per-sample
# single-cell files directly: many studies (e.g. GSE132771) ship only a
# `GSE..._RAW.tar`, with the loadable 10x triplets / h5 / h5ad living in each
# sample's own GSM suppl directory. When the series level yields nothing
# loadable, the manifest falls back to enumerating the series' samples (via
# getGEO()) and inventorying each GSM suppl directory. A GSM accession may also
# be passed directly to inventory or load a single sample.

# Classify a single supplementary filename into a (format, role) pair.
.classify_sc_file <- function(fname) {
    low <- tolower(fname)
    if (grepl("matrix\\.mtx(\\.gz)?$", low)) return(c(format = "10x_mtx", role = "matrix"))
    if (grepl("barcodes\\.tsv(\\.gz)?$", low)) return(c(format = "10x_mtx", role = "barcodes"))
    if (grepl("(features|genes)\\.tsv(\\.gz)?$", low)) return(c(format = "10x_mtx", role = "features"))
    if (grepl("\\.h5ad(\\.gz)?$", low)) return(c(format = "h5ad", role = "anndata"))
    if (grepl("\\.h5(\\.gz)?$", low)) return(c(format = "10x_h5", role = "matrix"))
    if (grepl("\\.loom(\\.gz)?$", low)) return(c(format = "loom", role = "matrix"))
    if (grepl("\\.(rds|rdata)(\\.gz)?$", low)) return(c(format = "rds", role = "object"))
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

# Formats readGEOSingleCell() can actually import: 10x Matrix Market, 10x HDF5,
# AnnData h5ad, and .rds holding a SingleCellExperiment or Seurat object
# (ADR-0006). loom is classified and reported but has no built-in reader (read
# it with its native package), so it is never "loadable".
.sc_readable <- c("10x_mtx", "10x_h5", "h5ad", "rds")

# The grouping key that defines a loadable unit. A 10x Matrix Market triplet is
# one unit per sample (its matrix/barcodes/features load together); every other
# format is one unit *per file*, because a study may ship several independent
# single-file objects (e.g. GSE161228's three whole-series .h5ad files) that must
# not be merged into a single bogus unit.
.sc_unit_key <- function(manifest) {
    ifelse(
        manifest$format == "10x_mtx",
        paste0("mtx|", manifest$sample),
        paste0("file|", manifest$fname)
    )
}

# Group a manifest into loadable units and assess completeness. Pure logic; used
# by geoSingleCellUnits().
.sc_units <- function(manifest) {
    cols <- c("unit", "sample", "platform", "format", "n_files", "status", "loadable")
    if (nrow(manifest) == 0) {
        out <- data.frame(matrix(nrow = 0, ncol = length(cols)))
        colnames(out) <- cols
        return(out)
    }
    has_plat <- "platform" %in% names(manifest)
    groups <- split(seq_len(nrow(manifest)), .sc_unit_key(manifest))
    parts <- Map(function(idx, k) {
        g <- manifest[idx, , drop = FALSE]
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
            unit = k,
            sample = g$sample[1],
            platform = if (has_plat) g$platform[1] else NA_character_,
            format = fmt, n_files = nrow(g), status = status,
            loadable = identical(status, "complete") && fmt %in% .sc_readable,
            stringsAsFactors = FALSE
        )
    }, groups, names(groups))
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out[order(out$sample, out$format), ]
}

#' Group a single-cell manifest into loadable units
#'
#' Collapses a \code{\link{geoSingleCellManifest}} into one row per unit and
#' reports completeness. A 10x Matrix Market unit groups a sample's matrix,
#' barcodes, and features files and is "complete" only when all three are
#' present; every other format is one unit per file. The \code{loadable} column
#' flags units a built-in reader can consume -- complete 10x Matrix Market, 10x
#' HDF5, and AnnData h5ad. loom and Seurat \code{.rds} are reported but not
#' loadable (read them with their native packages).
#'
#' @param manifest A data.frame returned by \code{geoSingleCellManifest()}.
#' @return A data.frame with columns \code{unit} (the grouping key),
#'   \code{sample}, \code{platform} (GPL, or NA), \code{format}, \code{n_files},
#'   \code{status}, and \code{loadable}.
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

# Inventory a single accession's suppl directory and classify the files. For a
# GSM, any filename that does not itself embed the GSM id still belongs to that
# sample, so fill NA sample ids with the accession.
.sc_manifest_for_accession <- function(GEO) {
    files <- getGEOSuppFiles(GEO, fetch_files = FALSE, quiet = TRUE)
    if (is.null(files) || nrow(files) == 0) {
        return(.classify_sc_files(character(0), character(0)))
    }
    m <- .classify_sc_files(files$fname, files$url)
    if (toupper(substr(GEO, 1, 3)) == "GSM") {
        m$sample[is.na(m$sample)] <- GEO
    }
    m
}

# Map the GSM accessions of a GSE to their platform (GPL). Uses getGEO() (Series
# Matrix path, GPL fetch skipped for speed): each returned object is one platform
# and its colData carries `platform_id` per sample. Returns a named character
# vector (names = GSM ids, values = GPL), empty on any failure. names() is also
# the GSM list used by the GSM-level manifest fallback.
.gse_gsm_platforms <- function(GEO) {
    empty <- stats::setNames(character(0), character(0))
    es <- tryCatch(
        getGEO(GEO, GSEMatrix = TRUE, getGPL = FALSE),
        error = function(e) NULL
    )
    if (is.null(es)) {
        return(empty)
    }
    if (!is.list(es)) {
        es <- list(es)
    }
    parts <- lapply(es, function(x) {
        ids <- tryCatch(colnames(x), error = function(e) NULL)
        if (is.null(ids)) {
            return(NULL)
        }
        plat <- tryCatch({
            cd <- SummarizedExperiment::colData(x)
            if ("platform_id" %in% colnames(cd)) {
                as.character(cd$platform_id)
            } else {
                rep(NA_character_, length(ids))
            }
        }, error = function(e) rep(NA_character_, length(ids)))
        stats::setNames(plat, ids)
    })
    # unname() so unlist() keeps the inner GSM-id names rather than prefixing
    # them with each list element's (file) name.
    out <- unlist(unname(parts))
    if (is.null(out)) {
        return(empty)
    }
    out <- out[!is.na(names(out)) & grepl("^GSM", names(out))]
    out[!duplicated(names(out))]
}

# Build a single-cell manifest from a set of GSM accessions by inventorying each
# sample's suppl directory and stacking the results.
.sc_manifest_from_gsms <- function(gsms) {
    parts <- lapply(gsms, function(g) {
        tryCatch(.sc_manifest_for_accession(g), error = function(e) NULL)
    })
    parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, parts)
    if (length(parts) == 0) {
        return(.classify_sc_files(character(0), character(0)))
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
}

#' Inventory the single-cell supplementary files of a GEO Series or Sample
#'
#' Lists the supplementary files attached to a GSE (or a single GSM) and
#' classifies each by single-cell format (10x Matrix Market triplet, 10x HDF5,
#' AnnData h5ad, loom, Seurat rds, tar archive, or other), extracting the GSM
#' sample id where present. This lets you see what a single-cell study contains
#' -- and how 10x triplets group by sample -- before downloading potentially
#' many gigabytes.
#'
#' For a GSE, the series-level suppl directory is inventoried first. Many
#' single-cell studies ship only a \code{GSE..._RAW.tar} there, with the
#' loadable per-sample files (10x triplets, h5, h5ad) living in each sample's
#' own GSM suppl directory. When the series level yields no loadable units, the
#' manifest falls back to enumerating the series' samples (via
#' \code{\link{getGEO}}) and inventorying each GSM suppl directory. Pass a GSM
#' accession to inventory just that one sample.
#'
#' No files are downloaded. The result feeds the single-cell readers (see
#' ADR-0004); reading itself uses Bioconductor importers (TENxIO, anndataR)
#' that are optional dependencies.
#'
#' The \code{platform} column (GPL accession per sample) is populated when the
#' manifest is built from the GSM level -- the common single-cell case, and the
#' one that matters, since a GSE can span multiple platforms (e.g. GSE132771
#' mixes mouse and human). It is the grouping used by
#' \code{getGEOSingleCell(by = "platform")}. It is \code{NA} for a single GSM,
#' and for whole-study files attached at the series level (which have no GSM).
#'
#' @param GEO A GEO Series (\code{"GSE..."}) or Sample (\code{"GSM..."})
#'   accession, e.g. "GSE132771" or "GSM3891612".
#' @param samples Optional character vector of GSM ids. For a GSE, restricts the
#'   inventory to these samples; when the series level has no loadable units this
#'   also avoids enumerating the whole series. Ignored when \code{GEO} is a GSM.
#' @return A data.frame with columns \code{fname}, \code{sample} (GSM id or NA),
#'   \code{platform} (GPL accession or NA), \code{format}, \code{role}, and
#'   \code{url}. Zero rows if nothing is found.
#' @seealso \code{\link{getGEOSuppFiles}}, \code{\link{getGEOSingleCell}}
#' @examples
#' \dontrun{
#'   geoSingleCellManifest("GSE132771")        # GSE: falls back to GSM level
#'   geoSingleCellManifest("GSM3891612")       # a single sample
#' }
#' @export
geoSingleCellManifest <- function(GEO, samples = NULL) {
    geotype <- toupper(substr(GEO, 1, 3))
    if (geotype == "GSM") {
        m <- .sc_manifest_for_accession(GEO)
        m$platform <- rep(NA_character_, nrow(m))
        return(m)
    }
    if (geotype != "GSE") {
        stop("geoSingleCellManifest() requires a GSE or GSM accession; got '",
            GEO, "'.", call. = FALSE)
    }
    # Series-level suppl files first. If they already hold loadable units (e.g.
    # whole-study .h5ad files), use them as-is; platform is not resolved here
    # (those files have no GSM, and resolving it would cost an extra download).
    series_m <- .sc_manifest_for_accession(GEO)
    if (any(.sc_units(series_m)$loadable %in% TRUE)) {
        if (!is.null(samples)) {
            series_m <- series_m[series_m$sample %in% samples, , drop = FALSE]
            rownames(series_m) <- NULL
        }
        series_m$platform <- rep(NA_character_, nrow(series_m))
        return(series_m)
    }
    # Otherwise fall back to per-sample (GSM) suppl directories. The GSM->GPL map
    # both enumerates the samples and supplies their platform.
    map <- .gse_gsm_platforms(GEO)
    gsms <- if (!is.null(samples)) samples else names(map)
    if (length(gsms) == 0) {
        series_m$platform <- rep(NA_character_, nrow(series_m))
        return(series_m)
    }
    gsm_m <- .sc_manifest_from_gsms(gsms)
    if (nrow(gsm_m) == 0) {
        series_m$platform <- rep(NA_character_, nrow(series_m))
        return(series_m)
    }
    gsm_m$platform <- unname(map[gsm_m$sample])
    rownames(gsm_m) <- NULL
    gsm_m
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

# HDF5-based readers need an actual (uncompressed) file. GEO often gzips them
# (e.g. GSE161228's *.h5ad.gz), so transparently decompress a .gz to a temp file
# and return that path; non-.gz paths pass through unchanged.
.maybe_gunzip <- function(path) {
    if (!grepl("\\.gz$", path)) {
        return(path)
    }
    dest <- tempfile(fileext = sub("^.*?(\\.[^.]+)\\.gz$", "\\1", basename(path)))
    R.utils::gunzip(path, destname = dest, remove = FALSE, overwrite = TRUE)
    dest
}

.read_sc_10x_mtx <- function(x) {
    .require_pkg("TENxIO", "10x Matrix Market data")
    BiocIO::import(TENxIO::TENxFileList(.arrange_10x(x)))
}
.read_sc_10x_h5 <- function(path) {
    .require_pkg("TENxIO", "10x HDF5 data")
    BiocIO::import(TENxIO::TENxH5(.maybe_gunzip(path)))
}
.read_sc_h5ad <- function(path) {
    .require_pkg("anndataR", "AnnData (h5ad) data")
    anndataR::read_h5ad(.maybe_gunzip(path), as = "SingleCellExperiment")
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

# When a (named) sample offers more than one loadable format, keep just one, by
# preference (richer/standard first). Units with no sample (NA) -- whole-study
# single files -- are each independent and pass through untouched (never collapse
# several NA-sample files into one).
.format_priority <- c("h5ad", "10x_h5", "10x_mtx")
.prefer_one_format <- function(load) {
    if (nrow(load) == 0) {
        return(load)
    }
    na_rows <- load[is.na(load$sample), , drop = FALSE]
    named <- load[!is.na(load$sample), , drop = FALSE]
    keep <- if (nrow(named) == 0) {
        named
    } else {
        do.call(rbind, lapply(split(named, named$sample), function(g) {
            if (nrow(g) == 1) {
                return(g)
            }
            ord <- order(match(g$format, .format_priority))
            g[ord[1], , drop = FALSE]
        }))
    }
    out <- rbind(keep, na_rows)
    rownames(out) <- NULL
    out[order(out$sample), , drop = FALSE]
}

# Pure unit-selection logic: split units into those to load and those skipped,
# honoring optional `samples` / `format` filters and one-format-per-sample.
# Units are identified by their `unit` key (so distinct whole-study files of the
# same format are tracked separately).
.select_sc_units <- function(units, samples = NULL, format = NULL) {
    load <- units[units$loadable %in% TRUE, , drop = FALSE]
    if (!is.null(samples)) {
        load <- load[load$sample %in% samples, , drop = FALSE]
    }
    if (!is.null(format)) {
        load <- load[load$format %in% format, , drop = FALSE]
    }
    load <- .prefer_one_format(load)
    skip <- units[!(units$unit %in% load$unit), , drop = FALSE]
    list(load = load, skip = skip)
}

# Download the files of one loadable unit (rows of the manifest, each carrying a
# full `url`) into destdir, skipping any already present, and return the local
# paths. Downloading by URL -- rather than re-listing a suppl directory -- works
# whether the files live at the series level or in a GSM suppl directory.
.download_sc_unit <- function(unit_files, destdir) {
    dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
    vapply(seq_len(nrow(unit_files)), function(i) {
        dest <- file.path(destdir, unit_files$fname[i])
        if (!file.exists(dest)) {
            downloadFile(unit_files$url[i], dest, quiet = TRUE)
        }
        dest
    }, character(1))
}

# cbind a list of SingleCellExperiments into one. Samples in a study often come
# from different references/platforms (e.g. GSE132771 mixes mouse and human) or
# CellRanger versions, so a naive cbind() fails two ways: the feature sets
# (rownames) differ, and even when they match the rowData *columns* differ (10x
# v2 genes.tsv -> ID,Symbol vs v3 features.tsv -> ID,Symbol,Type), which trips
# cbind's "rowData must be identical" check ("subscript contains invalid names").
# Align every object to the features common to all, in one consistent order, and
# give them one canonical rowData (the first sample's, restricted to the columns
# shared by all) so the rows match exactly before binding. Error clearly when
# there is no shared feature space to combine on.
.combine_sce <- function(results) {
    common <- Reduce(intersect, lapply(results, rownames))
    if (length(common) == 0) {
        stop(
            "Cannot combine: the samples share no common features (rownames). ",
            "They likely come from different platforms or genome references ",
            "(e.g. a study mixing organisms). Use `by = \"platform\"` to combine ",
            "within each platform, restrict with `samples=`, or `by = \"sample\"` ",
            "to keep them separate.",
            call. = FALSE
        )
    }
    dropped <- sum(vapply(results, function(x) nrow(x) - length(common), integer(1)))
    if (dropped > 0) {
        message(sprintf(
            "Combining on %d feature(s) shared by all %d samples (dropped %d non-shared feature row(s) across samples).",
            length(common), length(results), dropped
        ))
    }
    aligned <- lapply(results, function(x) x[common, ])
    # Impose one canonical rowData so heterogeneous per-sample feature annotation
    # (differing columns or values for the same feature) cannot block the cbind.
    common_cols <- Reduce(
        intersect,
        lapply(aligned, function(x) colnames(SummarizedExperiment::rowData(x)))
    )
    canon <- SummarizedExperiment::rowData(aligned[[1]])[, common_cols, drop = FALSE]
    aligned <- lapply(aligned, function(x) {
        SummarizedExperiment::rowData(x) <- canon
        x
    })
    do.call(SummarizedExperiment::cbind, aligned)
}

#' Download and read the single-cell data of a GEO Series or Sample
#'
#' High-level, best-effort convenience wrapper: inventories the GSE (or GSM)
#' (\code{\link{geoSingleCellManifest}}), groups files into loadable units
#' (\code{\link{geoSingleCellUnits}}), downloads each loadable unit, reads it
#' with \code{\link{readGEOSingleCell}}, and returns the results. It reports
#' which units it loads and which it skips.
#'
#' This handles common, well-structured layouts (clean per-sample 10x, h5ad, or
#' a saved object in \code{.rds}), including the very common case where the
#' series ships only a \code{_RAW.tar} and the per-sample files live in each GSM
#' suppl directory (the manifest falls back to the GSM level automatically). You
#' may also pass a single GSM accession to load just that sample. It does NOT
#' handle every GSE: loom files, files available \emph{only} inside a
#' \code{_RAW.tar} archive, and idiosyncratic layouts (e.g. a single combined
#' matrix for many samples) are out of scope -- use the manifest plus
#' \code{readGEOSingleCell()} directly for those.
#'
#' \strong{Grouping (\code{by}).} A GEO Series has two natural layers -- it can
#' span multiple platforms (GPLs), each holding many samples (GSMs) -- and the
#' platform is the feature-compatibility boundary (samples in one platform share
#' a feature space; across platforms they generally do not). \code{by} chooses
#' the return shape, and the shape is fixed by the argument (not the data):
#' \describe{
#'   \item{\code{"sample"} (default)}{a named list with one object per sample
#'     (or per whole-study file).}
#'   \item{\code{"platform"}}{a named list keyed by platform (GPL), each entry
#'     the samples of that platform combined into one object. The honest answer
#'     for a multi-platform study; a single-platform study yields a length-1
#'     list. Samples with unknown platform are returned individually.}
#'   \item{\code{"all"}}{a single object with every sample combined. Errors if
#'     the samples share no common features (e.g. a study mixing organisms) --
#'     use \code{"platform"} for those.}
#' }
#' Combining (for \code{"platform"}/\code{"all"}) restricts to the features
#' common to the group and reconciles per-sample feature annotation so binding
#' succeeds across CellRanger versions; whole-study single-file formats already
#' hold one object, so grouping is effectively a no-op for them.
#'
#' @param GEO A GEO Series (\code{"GSE..."}) or Sample (\code{"GSM..."})
#'   accession, e.g. "GSE132771" or "GSM3891612".
#' @param samples Optional character vector of GSM ids to restrict to. Ignored
#'   when \code{GEO} is itself a GSM.
#' @param format Optional format(s) to restrict to ("10x_mtx", "10x_h5",
#'   "h5ad", "rds").
#' @param by One of \code{"sample"} (default), \code{"platform"}, or
#'   \code{"all"} -- how to group the loaded samples into the return value. See
#'   Details.
#' @param as Output class, one of "SingleCellExperiment" (default) or "Seurat"
#'   (coerced at the boundary via the Seurat package, an optional dependency).
#' @param destdir Download destination directory.
#' @return Depends on \code{by}: a named list of objects per sample
#'   (\code{"sample"}); a named list of combined objects per platform
#'   (\code{"platform"}); or a single combined object (\code{"all"}). Each
#'   object is a \code{SingleCellExperiment}, or a \code{Seurat} object when
#'   \code{as = "Seurat"}.
#' @seealso \code{\link{geoSingleCellManifest}}, \code{\link{readGEOSingleCell}}
#' @examples
#' \dontrun{
#'   sce <- getGEOSingleCell("GSM3891612")                   # one sample
#'   per_sample <- getGEOSingleCell("GSE132771")             # list by GSM
#'   per_platform <- getGEOSingleCell("GSE132771", by = "platform")
#'   # -> list(GPL21103 = <mouse SCE>, GPL24676 = <human SCE>)
#' }
#' @export
getGEOSingleCell <- function(GEO, samples = NULL, format = NULL,
    by = c("sample", "platform", "all"),
    as = c("SingleCellExperiment", "Seurat"), destdir = tempdir()) {
    by <- match.arg(by)
    as <- match.arg(as)
    manifest <- geoSingleCellManifest(GEO, samples = samples)
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
    mkey <- .sc_unit_key(manifest)
    results <- list()
    for (i in seq_len(nrow(sel$load))) {
        u <- sel$load[i, ]
        unit_files <- manifest[mkey == u$unit, , drop = FALSE]
        label <- if (!is.na(u$sample)) u$sample else .sc_unit_label(unit_files$fname[1])
        message(sprintf("Loading %s (%s)...", label, u$format))
        local <- .download_sc_unit(unit_files, destdir)
        results[[label]] <- readGEOSingleCell(local, format = u$format)
    }
    if (by == "sample") {
        return(lapply(results, .as_sc_output, as = as))
    }
    if (by == "all") {
        combined <- if (length(results) == 1) results[[1]] else .combine_sce(results)
        return(.as_sc_output(combined, as))
    }
    # by == "platform": combine within each platform; samples of unknown platform
    # (NA) are kept individual (we cannot assert they share a feature space).
    plats <- sel$load$platform
    keys <- ifelse(is.na(plats), paste0("\r", names(results)), plats)
    groups <- split(seq_along(results), keys)
    out <- lapply(groups, function(idx) {
        grp <- results[idx]
        if (length(grp) == 1) grp[[1]] else .combine_sce(grp)
    })
    names(out) <- sub("^\r", "", names(out))
    lapply(out, .as_sc_output, as = as)
}

# A human-readable label for a whole-study (no-GSM) unit: the file name with the
# single-cell extension (and any .gz) stripped.
.sc_unit_label <- function(fname) {
    sub("\\.(mtx|tsv|h5ad|h5|rds|loom|tar)(\\.gz)?$", "", basename(fname))
}

# Persistent download cache backed by BiocFileCache (#171).
#
# Historically GEOquery "caches" by re-using files already present in `destdir`,
# which does not survive across sessions or guard partial downloads. This adds a
# proper persistent cache, keyed on the download URL. It is opt-in for one
# release cycle: set options(GEOquery.cache = TRUE) (or, per call, leave the
# default to preserve the historical destdir behavior). The cache location can
# be overridden with options(GEOquery.cache.path = ...).

.geo_cache_dir <- function() {
    getOption("GEOquery.cache.path", tools::R_user_dir("GEOquery", "cache"))
}

#' GEOquery download cache
#'
#' Return the \code{BiocFileCache} object that backs GEOquery's persistent
#' download cache (see \code{\link{clearGEOCache}}). The cache is used by the
#' download functions only when \code{options(GEOquery.cache = TRUE)} is set;
#' its location defaults to \code{tools::R_user_dir("GEOquery", "cache")} and can
#' be overridden with \code{options(GEOquery.cache.path = ...)}.
#'
#' @return A \code{BiocFileCache} object.
#' @seealso \code{\link{clearGEOCache}}
#' @export
geoCache <- function() {
    BiocFileCache::BiocFileCache(.geo_cache_dir(), ask = FALSE)
}

# Return the cached local path for `url`, or NULL if not cached.
.cache_lookup <- function(url, bfc = geoCache()) {
    q <- BiocFileCache::bfcquery(bfc, query = url, field = "rname", exact = TRUE)
    if (nrow(q) == 0) {
        return(NULL)
    }
    q$rpath[1]
}

# Add `path` to the cache under the key `url` (stored as a copy).
.cache_add <- function(url, path, bfc = geoCache()) {
    BiocFileCache::bfcadd(bfc, rname = url, fpath = path, rtype = "local", action = "copy")
    invisible(TRUE)
}

#' Clear the GEOquery download cache
#'
#' Remove all entries from the persistent download cache (see
#' \code{\link{geoCache}}).
#'
#' @return \code{NULL}, invisibly.
#' @seealso \code{\link{geoCache}}
#' @export
clearGEOCache <- function() {
    bfc <- geoCache()
    info <- BiocFileCache::bfcinfo(bfc)
    if (nrow(info) > 0) {
        BiocFileCache::bfcremove(bfc, info$rid)
    }
    invisible(NULL)
}

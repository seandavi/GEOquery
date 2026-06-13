# GEOquery download cache

Return the `BiocFileCache` object that backs GEOquery's persistent
download cache (see
[`clearGEOCache`](http://seandavi.github.io/GEOquery/reference/clearGEOCache.md)).
The cache is used by the download functions only when
`options(GEOquery.cache = TRUE)` is set; its location defaults to
`tools::R_user_dir("GEOquery", "cache")` and can be overridden with
`options(GEOquery.cache.path = ...)`.

## Usage

``` r
geoCache()
```

## Value

A `BiocFileCache` object.

## See also

[`clearGEOCache`](http://seandavi.github.io/GEOquery/reference/clearGEOCache.md)

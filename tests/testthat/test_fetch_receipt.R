# Synthetic bytes and mocked HTTP only. No recorded GEO responses.
receipt_fixture <- function(code) {
  directory <- tempfile("receipt-test-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  path <- file.path(directory, "receipt.json")
  old <- options(GEOquery.fetch.receipt = path, GEOquery.cache = FALSE)
  on.exit(options(old), add = TRUE)
  code(directory, path)
}

receipt_read <- function(path) {
  expect_true(file.exists(path))
  jsonlite::read_json(path, simplifyVector = FALSE)
}

receipt_gzip <- function(text) {
  path <- tempfile(fileext = ".gz")
  on.exit(unlink(path))
  connection <- gzfile(path, "wb")
  writeBin(charToRaw(paste0(paste(text, collapse = "\n"), "\n")), connection)
  close(connection)
  readBin(path, "raw", n = file.info(path)$size)
}

receipt_matrix <- function() {
  receipt_gzip(c(
    "!Series_title\tSynthetic series",
    "!Series_geo_accession\tGSE10",
    "!Series_summary\tSynthetic receipt fixture",
    "!Series_last_update_date\tJul 18 2026",
    "!Series_type\tExpression profiling by array",
    "!Sample_title\tA\tB",
    "!Sample_geo_accession\tGSM11\tGSM12",
    "!Sample_platform_id\tGPL10\tGPL10",
    "!series_matrix_table_begin",
    "\"ID_REF\"\t\"GSM11\"\t\"GSM12\"",
    "\"probe_1\"\t1.5\t2.5",
    "\"probe_2\"\t3.5\t4.5",
    "!series_matrix_table_end"
  ))
}

receipt_gpl <- function() {
  receipt_gzip(c(
    "^PLATFORM = GPL10",
    "!Platform_title = Synthetic platform",
    "!Platform_geo_accession = GPL10",
    "!Platform_last_update_date = Jul 17 2026",
    "#ID = Probe ID",
    "#GENE = Gene label",
    "!platform_table_begin",
    "ID\tGENE",
    "probe_1\tA",
    "probe_2\tB",
    "!platform_table_end"
  ))
}

receipt_index <- function(names) {
  charToRaw(paste0("<a href=\"", names, "\">", names, "</a>", collapse = ""))
}

receipt_response <- function(bytes, path = NULL) {
  # httr2 mocked responses bypass req_perform(path=)'s streaming writer.
  # Materialize exactly the authored response bytes for download fixtures.
  if (!is.null(path)) {
    writeBin(bytes, path)
  }
  httr2::response(
    200, body = bytes,
    headers = list(`last-modified` = "Sat, 18 Jul 2026 12:00:00 GMT")
  )
}

test_that("getGEO receipts include every matrix and nested GPL only once", {
  receipt_fixture(function(directory, path) {
    files <- c(
      "GSE10-GPL10_series_matrix.txt.gz", "GSE10-GPL11_series_matrix.txt.gz"
    )
    requests <- character()
    matrix_bytes <- receipt_matrix()
    gpl_bytes <- receipt_gpl()
    mock <- function(req) {
      requests <<- c(requests, req$url)
      if (grepl("/matrix/$", req$url)) {
        receipt_response(receipt_index(files))
      } else if (grepl("acc=GPL10", req$url)) {
        receipt_response(gpl_bytes, file.path(directory, "GPL10.soft.gz"))
      } else {
        receipt_response(matrix_bytes, file.path(directory, basename(req$url)))
      }
    }
    value <- httr2::with_mocked_responses(
      mock, getGEO(
        "GSE10", destdir = directory, returnType = "SummarizedExperiment"
      )
    )
    expect_length(requests, 4L)
    expect_length(value, 2L)
    expect_true(all(vapply(
      value, methods::is, logical(1), "SummarizedExperiment"
    )))
    document <- receipt_read(path)
    expect_equal(document$schema_version, 1L)
    expect_identical(document$operation, "getGEO")
    expect_length(document$resources, 3L)
    expect_equal(
      vapply(document$resources, `[[`, "", "accession"),
      c("GSE10", "GPL10", "GSE10")
    )
    dates <- vapply(document$resources, `[[`, "", "geo_last_update")
    expect_equal(dates, c("Jul 18 2026", "Jul 17 2026", "Jul 18 2026"))
    for (resource in document$resources) {
      expect_named(resource, c(
        "accession", "url", "endpoint", "expected_url", "file", "source",
        "bytes", "md5", "http_last_modified", "geo_last_update",
        "geoquery_version", "timestamp"
      ))
      file <- file.path(directory, resource$file)
      expect_identical(resource$source, "download")
      expect_equal(resource$bytes, unname(file.info(file)$size))
      expect_identical(resource$md5, unname(tools::md5sum(file)))
      expect_identical(
        resource$http_last_modified, "Sat, 18 Jul 2026 12:00:00 GMT"
      )
      expect_identical(
        resource$geoquery_version, as.character(packageVersion("GEOquery"))
      )
      expect_match(resource$timestamp, "^20[0-9-]+T[0-9:]+Z$")
      expect_null(resource$expected_url)
    }
    expect_null(GEOquery:::.geo_receipt_state$active)
  })
})

test_that("existing matrices record local bytes and unknown HTTP metadata", {
  receipt_fixture(function(directory, path) {
    filename <- "GSE10_series_matrix.txt.gz"
    writeBin(receipt_matrix(), file.path(directory, filename))
    requests <- 0L
    mock <- function(req) {
      requests <<- requests + 1L
      expect_match(req$url, "/matrix/$")
      receipt_response(receipt_index(filename))
    }
    httr2::with_mocked_responses(
      mock, getGEO("GSE10", destdir = directory, getGPL = FALSE)
    )
    expect_equal(requests, 1L)
    resource <- receipt_read(path)$resources[[1L]]
    expect_identical(resource$source, "local")
    expect_null(resource$url)
    expect_null(resource$endpoint)
    expect_null(resource$http_last_modified)
    expect_match(resource$expected_url, "GSE10_series_matrix.txt.gz$")
    expect_identical(resource$geo_last_update, "Jul 18 2026")
  })
})

test_that("persistent cache hits record their key and no invented HTTP date", {
  receipt_fixture(function(directory, path) {
    filename <- "GSE10_series_matrix.txt.gz"
    url <- paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE10/matrix/", filename
    )
    fixture <- file.path(directory, "fixture.gz")
    writeBin(receipt_matrix(), fixture)
    old <- options(
      GEOquery.cache = TRUE, GEOquery.cache.path = file.path(directory, "cache")
    )
    on.exit(options(old), add = TRUE)
    GEOquery:::.cache_add(url, fixture)
    requests <- 0L
    mock <- function(req) {
      requests <<- requests + 1L
      expect_match(req$url, "/matrix/$")
      receipt_response(receipt_index(filename))
    }
    httr2::with_mocked_responses(
      mock, getGEO("GSE10", destdir = directory, getGPL = FALSE)
    )
    expect_equal(requests, 1L)
    resource <- receipt_read(path)$resources[[1L]]
    expect_identical(resource$source, "cache")
    expect_identical(resource$url, url)
    expect_null(resource$http_last_modified)
    expect_identical(resource$md5, unname(tools::md5sum(fixture)))
  })
})

test_that("local inputs retain their return type and parsed metadata", {
  receipt_fixture(function(directory, path) {
    fixture <- file.path(directory, "local.gz")
    writeBin(receipt_matrix(), fixture)
    value <- httr2::with_mocked_responses(
      function(req) stop("Unexpected network access"),
      getGEO(filename = fixture, getGPL = FALSE, returnType = "ExpressionSet")
    )
    expect_s4_class(value, "ExpressionSet")
    expect_equal(unname(Biobase::exprs(value)[1, 1]), 1.5)
    resource <- receipt_read(path)$resources[[1L]]
    expect_identical(resource$accession, "GSE10")
    expect_identical(resource$geo_last_update, "Jul 18 2026")
    expect_null(resource$url)
    expect_null(resource$expected_url)
  })
})

test_that("supplementary receipts cover filtered downloads and local reuse", {
  receipt_fixture(function(directory, path) {
    files <- c("GSE10_A.txt", "GSE10_B.txt", "GSE10_RAW.tar")
    writeLines("existing", file.path(directory, files[[2L]]))
    requests <- character()
    mock <- function(req) {
      requests <<- c(requests, req$url)
      if (grepl("/suppl/$", req$url)) {
        receipt_response(receipt_index(files))
      } else {
        receipt_response(
          charToRaw("new bytes"), file.path(directory, files[[1L]])
        )
      }
    }
    value <- httr2::with_mocked_responses(mock, getGEOSuppFiles(
      "GSE10", makeDirectory = FALSE, baseDir = directory,
      filter_regex = "\\.txt$", quiet = TRUE
    ))
    expect_equal(value$fname, files[1:2])
    expect_length(requests, 2L)
    document <- receipt_read(path)
    expect_identical(document$operation, "getGEOSuppFiles")
    expect_length(document$resources, 2L)
    expect_identical(document$resources[[1L]]$source, "download")
    expect_identical(document$resources[[2L]]$source, "local")
    for (resource in document$resources) {
      expect_identical(resource$accession, "GSE10")
      expect_null(resource$geo_last_update)
      expect_identical(
        resource$md5,
        unname(tools::md5sum(file.path(directory, resource$file)))
      )
    }
  })
})

test_that("listing supplementary files leaves receipts disabled", {
  receipt_fixture(function(directory, path) {
    value <- httr2::with_mocked_responses(
      function(req) receipt_response(receipt_index("GSE10_A.txt")),
      getGEOSuppFiles("GSE10", baseDir = directory, fetch_files = FALSE)
    )
    expect_equal(value$fname, "GSE10_A.txt")
    expect_false(file.exists(path))
    expect_null(GEOquery:::.geo_receipt_state$active)
  })
})

test_that("SOFT receipts redact reviewer tokens and retain parsed GEO dates", {
  receipt_fixture(function(directory, path) {
    body <- paste(c(
      "^SAMPLE = GSM10", "!Sample_title = Synthetic sample",
      "!Sample_geo_accession = GSM10", "!Sample_last_update_date = Jul 16 2026",
      "#ID_REF = Probe", "#VALUE = Measurement", "!sample_table_begin",
      "ID_REF\tVALUE", "probe_1\t2", "!sample_table_end"
    ), collapse = "\n")
    seen <- NULL
    value <- httr2::with_mocked_responses(function(req) {
      seen <<- req$url
      receipt_response(
        charToRaw(paste0(body, "\n")), file.path(directory, "GSM10.soft")
      )
    }, getGEO("GSM10", token = "SYNTHETIC_REVIEWER_VALUE", destdir = directory))
    expect_s4_class(value, "GSM")
    expect_match(seen, "token=SYNTHETIC_REVIEWER_VALUE")
    resource <- receipt_read(path)$resources[[1L]]
    expect_identical(resource$geo_last_update, "Jul 16 2026")
    expect_match(resource$url, "acc=GSM10")
    expect_false(grepl(
      "token|SYNTHETIC_REVIEWER_VALUE", paste(readLines(path), collapse = "")
    ))
  })
})

test_that("URL sanitizing drops credentials, fragments and unknown fields", {
  url <- paste0(
    "https://user:fake@example.org/file?acc=GSE10",
    "&token=FAKE&X-Amz-Signature=FAKE#FAKE"
  )
  expect_identical(
    GEOquery:::.receipt_url(url), "https://example.org/file?acc=GSE10"
  )
})

test_that("accession extraction skips numbered FTP directory stubs", {
  urls <- c(
    paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",
           "GSE123nnn/GSE123456/matrix/GSE123456.txt.gz"),
    paste0("https://ftp.ncbi.nlm.nih.gov/geo/platforms/",
           "GPL22nnn/GPL22222/annot/GPL22222.annot.gz"),
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM12345&form=text"
  )
  expect_equal(
    vapply(urls, GEOquery:::.receipt_accession, "", USE.NAMES = FALSE),
    c("GSE123456", "GPL22222", "GSM12345")
  )
  expect_null(GEOquery:::.receipt_accession("https://example.org/GSE123nnn/"))
})

test_that("receipts are disabled by default and do not change parsed values", {
  receipt_fixture(function(directory, path) {
    fixture <- file.path(directory, "matrix.gz")
    writeBin(receipt_matrix(), fixture)
    options(GEOquery.fetch.receipt = NULL)
    original <- getGEO(
      filename = fixture, getGPL = FALSE, returnType = "ExpressionSet"
    )
    expect_false(file.exists(path))
    options(GEOquery.fetch.receipt = path)
    recorded <- getGEO(
      filename = fixture, getGPL = FALSE, returnType = "ExpressionSet"
    )
    expect_equal(recorded, original)
    expect_length(receipt_read(path)$resources, 1L)
  })
})

test_that("invalid or occupied destinations fail and preserve files", {
  receipt_fixture(function(directory, path) {
    invalid_paths <- list(
      TRUE, 1, NA_character_, "", c("a", "b"),
      file.path(directory, "absent", "r.json")
    )
    for (invalid in invalid_paths) {
      options(GEOquery.fetch.receipt = invalid)
      expect_error(getGEO("GSM10"), class = "geoquery_receipt_error")
      expect_null(GEOquery:::.geo_receipt_state$active)
    }
    writeLines("preserve me", path)
    options(GEOquery.fetch.receipt = path)
    expect_error(
      getGEO("GSM10"), "already exists", class = "geoquery_receipt_error"
    )
    expect_identical(readLines(path), "preserve me")
  })
})

test_that("a dangling symlink is an occupied receipt destination", {
  receipt_fixture(function(directory, path) {
    target <- file.path(directory, "absent.json")
    if (!isTRUE(suppressWarnings(file.symlink(target, path)))) {
      skip("This filesystem does not support creating symlinks")
    }
    expect_error(
      getGEO("GSM10"), "already exists", class = "geoquery_receipt_error"
    )
    expect_identical(Sys.readlink(path), target)
    expect_false(file.exists(target))
    expect_null(GEOquery:::.geo_receipt_state$active)
  })
})

test_that("HTTP and parser errors leave no receipt and clear the collector", {
  receipt_fixture(function(directory, path) {
    httr2::with_mocked_responses(
      function(req) httr2::response(404),
      expect_error(
        getGEO("GSM10", destdir = directory), class = "geoquery_download_error"
      )
    )
    expect_false(file.exists(path))
    expect_null(GEOquery:::.geo_receipt_state$active)
    bad <- file.path(directory, "invalid.soft")
    writeLines(c(
      "!Series_title\tmalformed", "!series_matrix_table_begin",
      "!series_matrix_table_begin"
    ), bad)
    expect_error(getGEO(filename = bad), "parsing failed")
    expect_false(file.exists(path))
    expect_null(GEOquery:::.geo_receipt_state$active)
    writeBin(receipt_matrix(), bad)
    expect_s4_class(
      getGEO(filename = bad, getGPL = FALSE), "SummarizedExperiment"
    )
    expect_length(receipt_read(path)$resources, 1L)
  })
})

test_that("a receipt preserves downloaded data at a colliding output path", {
  receipt_fixture(function(directory, path) {
    collision <- file.path(directory, "GSE10_A.txt")
    options(GEOquery.fetch.receipt = collision)
    httr2::with_mocked_responses(function(req) {
      if (grepl("/suppl/$", req$url)) {
        receipt_response(receipt_index("GSE10_A.txt"))
      } else {
        receipt_response(charToRaw("preserve downloaded bytes"), collision)
      }
    }, expect_error(
      getGEOSuppFiles("GSE10", baseDir = directory, makeDirectory = FALSE),
      "already exists", class = "geoquery_receipt_error"
    ))
    expect_identical(
      readChar(collision, file.info(collision)$size),
      "preserve downloaded bytes"
    )
    expect_null(GEOquery:::.geo_receipt_state$active)
    expect_false(any(grepl(
      "^\\.geoquery-receipt-", list.files(directory, all.files = TRUE)
    )))
  })
})

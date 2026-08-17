# Offline tests for the shared HTTP layer (#173). httr2's own mocking is used,
# so no network and no extra dependency.

test_that(".geo_request retries transport failures, not just bad statuses", {
    # NCBI refuses connections under load; httr2 does not retry those unless
    # retry_on_failure is set. Mocked errors bypass httr2's retry loop, so the
    # policy is asserted directly rather than exercised.
    pol <- GEOquery:::.geo_request("https://example.com/x")$policies
    expect_true(pol$retry_on_failure)
    expect_equal(pol$retry_max_tries, 3)
})

test_that("getDirListing parses GEO-prefixed hrefs from a mocked index (#173)", {
    html <- paste0(
        "<html><body>",
        "<a href=\"GSE1_RAW.tar\">GSE1_RAW.tar</a>",
        "<a href=\"filelist.txt\">filelist.txt</a>",
        "</body></html>"
    )
    httr2::with_mocked_responses(
        function(req) httr2::response(200, body = charToRaw(html)),
        {
            fnames <- GEOquery:::getDirListing("https://ftp.ncbi.nlm.nih.gov/x/")
            expect_true("GSE1_RAW.tar" %in% fnames)
            # filtered: only names starting with "G" are returned
            expect_false("filelist.txt" %in% fnames)
        }
    )
})

test_that("downloadFile succeeds on a mocked 200 response (#173)", {
    res <- httr2::with_mocked_responses(
        function(req) httr2::response(200, body = charToRaw("payload")),
        GEOquery:::downloadFile("https://example.com/f", tempfile())
    )
    expect_true(res)
})

test_that("downloadFile raises a typed geoquery_download_error on HTTP failure (#173, #170)", {
    httr2::with_mocked_responses(
        function(req) httr2::response(404),
        expect_error(
            GEOquery:::downloadFile("https://example.com/missing", tempfile()),
            class = "geoquery_download_error"
        )
    )
})

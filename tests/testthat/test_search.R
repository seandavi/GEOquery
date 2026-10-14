# searchFieldsGEO() wraps rentrez::entrez_db_searchable("gds"), an NCBI EInfo
# call. The offline test pins the reshaping against a canned EInfo result; the
# live test is gated like the rest of the network suite, because NCBI sometimes
# answers EInfo with an HTML error page, which failed R CMD check on the
# Bioconductor builders (#169 strategy).

fake_einfo_entry <- function(name, full_name, description) {
  structure(
    list(
      Name = name, FullName = full_name, Description = description,
      TermCount = "0", IsDate = "N", IsNumerical = "N",
      SingleToken = "N", Hierarchy = "N", IsHidden = "N"
    ),
    class = c("eInfoEntry", "list")
  )
}

test_that("searchFieldsGEO() reshapes EInfo fields into one row per field", {
  fake <- structure(
    list(
      ALL = fake_einfo_entry("ALL", "All Fields", "All terms"),
      ORGN = fake_einfo_entry("ORGN", "Organism", "Scientific name")
    ),
    class = c("eInfoSearch", "eInfoList", "list")
  )
  local_mocked_bindings(
    entrez_db_searchable = function(db, ...) {
      expect_identical(db, "gds")
      fake
    },
    .package = "rentrez"
  )

  result <- searchFieldsGEO()

  expect_s3_class(result, "data.frame")
  expect_identical(nrow(result), 2L)
  expect_identical(rownames(result), c("1", "2"))
  expect_true(
    all(c("Name", "Description", "FullName") %in% colnames(result))
  )
  expect_identical(unname(unlist(result$Name)), c("ALL", "ORGN"))
})

test_that("returns_correct_structure_when_geo_accessible", {
  skip_if_no_integration()
  result <- searchFieldsGEO()

  expect_s3_class(result, "data.frame")
  expect_true(
    all(c("Name", "Description", "FullName") %in%
      colnames(result))
  )
})

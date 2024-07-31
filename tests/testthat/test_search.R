test_that("returns_correct_structure_when_geo_accessible", {
  result <- searchFieldsGEO()

  expect_s3_class(result, "data.frame")
  expect_true(all(c("Field", "Description", "FullName") %in%
    colnames(result)))
})

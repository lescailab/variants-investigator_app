# Load application support files into testing environment
# File: tests/testthat/test-app-function.R
library(shinytest2)

test_that("variantsinvestigator sql initial values are consistent", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()

  app <- AppDriver$new(vcfToSQL(), name = "variantsinvestigator-vcftosql")

  expect_setequal(names(app$get_values()$input), c("run_convert", "sqlite_file", "vcf_file"))
})

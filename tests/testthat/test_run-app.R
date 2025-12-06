# Load application support files into testing environment
# File: tests/testthat/test-app-function.R
library(shinytest2)

test_that("variantsinvestigator basic app test", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()

  app <- AppDriver$new(VariantsInvestigator(), name = "variantsinvestigator-app")

  # We check for the presence of key inputs instead of a full snapshot
  # because some widgets (like pickerInput) can be flaky across environments
  values <- app$get_values()
  expect_true("impact" %in% names(values$input))
  expect_true("dark_mode" %in% names(values$input))
  expect_true("gnomad" %in% names(values$input))
})

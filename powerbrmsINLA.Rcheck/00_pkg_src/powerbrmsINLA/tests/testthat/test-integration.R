# tests/testthat/test-integration.R
# Integration tests

library(testthat)
library(brms)
library(ggplot2)
library(dplyr)

test_that("Package integration works end-to-end", {
  # Test that we can create a simple data generator
  generator <- .auto_data_generator(
    formula = y ~ treatment,
    effect_name = "treatment",
    family = gaussian()
  )

  expect_true(is.function(generator))

  # Test data generation
  data_test <- generator(n = 20, effect = c(treatment = 0.5))
  expect_s3_class(data_test, "data.frame")
  expect_equal(nrow(data_test), 20)
  expect_true("y" %in% names(data_test))
  expect_true("treatment" %in% names(data_test))
})

test_that("Formula conversion pipeline works", {
  # Test complete pipeline from brms formula to INLA
  original_formula <- outcome ~ treatment + age + (1 + time | subject)

  # Parse random effects
  re_specs <- .parse_re_terms(original_formula)
  expect_length(re_specs, 1)
  expect_equal(re_specs[[1]]$group, "subject")
  expect_equal(re_specs[[1]]$slope, "time")

  # Convert to INLA formula
  conversion <- .brms_to_inla_formula2(original_formula)
  expect_s3_class(conversion$inla_formula, "formula")
  expect_length(conversion$re_specs, 1)

  # Check that INLA formula contains expected components
  formula_text <- as.character(conversion$inla_formula)
  expect_true(any(grepl("treatment", formula_text)))
  expect_true(any(grepl("age", formula_text)))
  expect_true(any(grepl("f\\(", formula_text)))  # Random effects
})

test_that("Error handling works consistently", {
  # Test that errors are handled gracefully across functions

  # Invalid prior specification
  invalid_priors <- data.frame(
    class = "b",
    coef = "treatment",
    prior = "invalid_distribution(0, 1)",
    stringsAsFactors = FALSE
  )
  result <- .map_brms_priors_to_inla(invalid_priors)
  expect_type(result, "list")  # Should not error, just ignore invalid priors

  # Test that the function properly validates input length
  expect_error(
    beta_weights_on_grid(c(0.5), mode = 0.5, n = 5),
    "length\\(effects\\) >= 2"
  )
})

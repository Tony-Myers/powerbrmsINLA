# tests/testthat/test-utils.R

library(testthat)
library(brms)

test_that("Coalesce operator works correctly", {
  expect_equal(NULL %||% "default", "default")
  expect_equal("value" %||% "default", "value")
  expect_equal(0 %||% "default", 0)
  expect_equal(FALSE %||% "default", FALSE)
})

test_that("ggplot2 line argument detection works", {
  arg <- .gg_line_arg()
  expect_true(arg %in% c("size", "linewidth"))
  expect_type(arg, "character")
  expect_length(arg, 1)
})

test_that("Family mapping works correctly", {
  # Test gaussian family
  result <- .to_inla_family(gaussian())
  expect_equal(result$inla, "gaussian")
  expect_equal(result$brms, "gaussian")

  # Test binomial family
  result <- .to_inla_family(binomial())
  expect_equal(result$inla, "binomial")
  expect_equal(result$brms, "binomial")

  # Test character input - create mock family object instead of using string directly
  result <- .to_inla_family("poisson")
  expect_equal(result$inla, "poisson")
  expect_equal(result$brms, "poisson")

  # Test student_t mapping
  result <- .to_inla_family("student")
  expect_equal(result$inla, "T")
})

test_that("Random effects parsing works", {
  # Simple random intercept
  formula1 <- y ~ x + (1 | group)
  re_specs <- .parse_re_terms(formula1)
  expect_length(re_specs, 1)
  expect_equal(re_specs[[1]]$group, "group")
  expect_true(re_specs[[1]]$has_intercept)
  expect_null(re_specs[[1]]$slope)

  # Random slope
  formula2 <- y ~ x + (1 + time | subject)
  re_specs <- .parse_re_terms(formula2)
  expect_length(re_specs, 1)
  expect_equal(re_specs[[1]]$group, "subject")
  expect_true(re_specs[[1]]$has_intercept)
  expect_equal(re_specs[[1]]$slope, "time")

  # No random effects
  formula3 <- y ~ x + z
  re_specs <- .parse_re_terms(formula3)
  expect_length(re_specs, 0)
})

test_that("brms to INLA formula conversion works", {
  # Simple formula
  formula1 <- y ~ x + (1 | group)
  result <- .brms_to_inla_formula2(formula1)
  expect_s3_class(result$inla_formula, "formula")
  expect_length(result$re_specs, 1)

  # Formula with multiple predictors
  formula2 <- outcome ~ treatment + age + (1 + time | subject)
  result <- .brms_to_inla_formula2(formula2)
  expect_s3_class(result$inla_formula, "formula")
  expect_length(result$re_specs, 1)
  expect_equal(result$re_specs[[1]]$slope, "time")

  # Check that formula text contains expected elements
  formula_text <- as.character(result$inla_formula)
  expect_true(any(grepl("treatment", formula_text)))
  expect_true(any(grepl("age", formula_text)))
  expect_true(any(grepl("f\\(", formula_text)))
})

test_that("Prior mapping works correctly", {
  # Test with NULL priors
  result <- .map_brms_priors_to_inla(NULL)
  expect_type(result, "list")
  expect_named(result, c("control_fixed", "hyper_by_re"))

  # Test with normal prior
  normal_priors <- data.frame(
    class = c("Intercept", "b"),
    coef = c("", "treatment"),
    prior = c("normal(0, 1)", "normal(0, 0.5)"),
    stringsAsFactors = FALSE
  )
  result <- .map_brms_priors_to_inla(normal_priors)
  expect_type(result$control_fixed, "list")
  expect_equal(result$control_fixed$mean.intercept, 0)
  expect_equal(result$control_fixed$prec.intercept, 1)
  expect_equal(result$control_fixed$mean$treatment, 0)
  expect_equal(result$control_fixed$prec$treatment, 4)
})

test_that("Wilson confidence interval stopping rule works", {
  # Test basic functionality
  result <- .should_stop_binom(8, 10, target = 0.8, margin = 0.1)
  expect_type(result, "list")
  expect_named(result, c("stop", "ci"))
  expect_type(result$stop, "logical")
  expect_length(result$ci, 2)
  expect_true(result$ci[1] <= result$ci[2])

  # Test edge cases
  result_zero <- .should_stop_binom(0, 0, target = 0.8)
  expect_false(result_zero$stop)
  expect_equal(result_zero$ci, c(0, 1))

  # Test clear stopping condition
  result_stop <- .should_stop_binom(10, 10, target = 0.5, margin = 0.1)
  expect_true(result_stop$stop)

  result_stop_low <- .should_stop_binom(0, 10, target = 0.5, margin = 0.1)
  expect_true(result_stop_low$stop)
})

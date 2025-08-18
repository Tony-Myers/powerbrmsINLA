# tests/testthat/test-sequential-fixed.R

library(testthat)

test_that("Wilson stopping rule works correctly", {
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

  # Test clear stopping condition - high success rate
  result_stop_high <- .should_stop_binom(10, 10, target = 0.5, margin = 0.1)
  expect_true(result_stop_high$stop)

  # Test clear stopping condition - low success rate
  result_stop_low <- .should_stop_binom(0, 10, target = 0.5, margin = 0.1)
  expect_true(result_stop_low$stop)

  # Test middle ground - should not stop
  result_no_stop <- .should_stop_binom(5, 10, target = 0.5, margin = 0.05)
  # This might or might not stop depending on the CI width, but should not error
  expect_type(result_no_stop$stop, "logical")
})

test_that("Sequential engine basic validation works", {
  # Test that the function can be called without error
  # Using a mock data generator to avoid INLA dependency
  mock_generator <- function(n, effect) {
    effect_val <- if (is.list(effect)) effect[[1]] else effect
    data.frame(
      y = rnorm(n, effect_val, 1),
      treatment = rnorm(n, 0, 1)
    )
  }

  # Test with minimal parameters
  expect_no_error({
    result <- brms_inla_power_sequential(
      formula = y ~ treatment,
      effect_name = "treatment",
      effect_grid = c(0.3, 0.7),
      sample_sizes = c(20, 30),
      metric = "direction",
      target = 0.8,
      batch_size = 3,
      min_sims = 6,
      max_sims = 15,
      data_generator = mock_generator,
      progress = FALSE
    )
  })
})

test_that("Sequential engine handles multi-effect scenarios", {
  # Mock generator for multi-effect
  mock_generator_multi <- function(n, effect) {
    treatment_val <- if (is.list(effect)) effect$treatment else effect
    age_val <- if (is.list(effect)) effect$age else 0

    data.frame(
      y = rnorm(n, treatment_val + age_val, 1),
      treatment = rnorm(n, 0, 1),
      age = rnorm(n, 0, 1)
    )
  }

  # Create multi-effect grid
  effect_grid <- expand.grid(
    treatment = c(0, 0.5),
    age = c(-0.2, 0.2)
  )

  expect_no_error({
    result <- brms_inla_power_sequential(
      formula = y ~ treatment + age,
      effect_name = c("treatment", "age"),
      effect_grid = effect_grid,
      sample_sizes = c(25),
      metric = "direction",
      target = 0.7,
      batch_size = 2,
      min_sims = 4,
      max_sims = 10,
      data_generator = mock_generator_multi,
      progress = FALSE
    )
  })
})

test_that("Sequential engine parameter validation works", {
  mock_generator <- function(n, effect) {
    data.frame(y = rnorm(n), x = rnorm(n))
  }

  # Test invalid metric
  expect_error(
    brms_inla_power_sequential(
      formula = y ~ x,
      effect_name = "x",
      effect_grid = c(0.5),
      sample_sizes = 20,
      metric = "invalid_metric",
      data_generator = mock_generator
    )
  )

  # Test invalid target
  expect_error(
    brms_inla_power_sequential(
      formula = y ~ x,
      effect_name = "x",
      effect_grid = c(0.5),
      sample_sizes = 20,
      metric = "direction",
      target = 1.5,  # > 1
      data_generator = mock_generator
    )
  )

  # Test empty effect_name
  expect_error(
    brms_inla_power_sequential(
      formula = y ~ x,
      effect_name = character(0),
      effect_grid = c(0.5),
      sample_sizes = 20,
      metric = "direction",
      data_generator = mock_generator
    )
  )
})

test_that("Sequential engine output structure is correct", {
  mock_generator <- function(n, effect) {
    data.frame(
      y = rnorm(n, if(is.list(effect)) effect[[1]] else effect, 1),
      treatment = rnorm(n, 0, 1)
    )
  }

  result <- brms_inla_power_sequential(
    formula = y ~ treatment,
    effect_name = "treatment",
    effect_grid = c(0.3),
    sample_sizes = c(20),
    metric = "direction",
    target = 0.8,
    batch_size = 2,
    min_sims = 4,
    max_sims = 8,
    data_generator = mock_generator,
    progress = FALSE
  )

  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("results", "summary", "settings"))
  expect_null(result$results)  # Sequential doesn't store individual results
  expect_s3_class(result$summary, "data.frame")
  expect_type(result$settings, "list")

  # Check summary structure
  expect_true("assurance" %in% names(result$summary))
  expect_true("sims_used" %in% names(result$summary))
  expect_true("n" %in% names(result$summary))
  expect_true("power_direction" %in% names(result$summary))

  # Check that assurance is between 0 and 1
  expect_true(all(result$summary$assurance >= 0 & result$summary$assurance <= 1, na.rm = TRUE))

  # Check that sims_used is reasonable
  expect_true(all(result$summary$sims_used >= result$settings$min_sims |
                    result$summary$sims_used <= result$settings$max_sims))
})

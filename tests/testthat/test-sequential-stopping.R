# tests/testthat/test-sequential-stopping.R
# Tests for the Wilson stopping rule and brms_inla_power_sequential().
# (Renamed from sequencial-test.R, which did not match testthat's test-*
#  pattern and was never executed.)

library(testthat)

test_that("Wilson stopping rule works correctly", {
  result <- .should_stop_binom(8, 10, target = 0.8, margin = 0.1)
  expect_type(result, "list")
  expect_named(result, c("stop", "ci"))
  expect_type(result$stop, "logical")
  expect_length(result$ci, 2)
  expect_true(result$ci[1] <= result$ci[2])

  # Edge case: no trials yet
  result_zero <- .should_stop_binom(0, 0, target = 0.8)
  expect_false(result_zero$stop)
  expect_equal(result_zero$ci, c(0, 1))

  # Clear stopping condition - high success rate
  result_stop_high <- .should_stop_binom(10, 10, target = 0.5, margin = 0.1)
  expect_true(result_stop_high$stop)

  # Clear stopping condition - low success rate
  result_stop_low <- .should_stop_binom(0, 10, target = 0.5, margin = 0.1)
  expect_true(result_stop_low$stop)

  # Middle ground - should not error
  result_no_stop <- .should_stop_binom(5, 10, target = 0.5, margin = 0.05)
  expect_type(result_no_stop$stop, "logical")
})

test_that("Sequential engine parameter validation works without INLA", {
  mock_generator <- function(n, effect) {
    data.frame(y = rnorm(n), x = rnorm(n))
  }

  # Invalid metric (errors at match.arg, before the INLA dependency check)
  expect_error(
    brms_inla_power_sequential(
      formula = y ~ x, effect_name = "x",
      effect_grid = c(0.5), sample_sizes = 20,
      metric = "invalid_metric", data_generator = mock_generator
    )
  )

  # Invalid target (must lie strictly in (0, 1))
  expect_error(
    brms_inla_power_sequential(
      formula = y ~ x, effect_name = "x",
      effect_grid = c(0.5), sample_sizes = 20,
      metric = "direction", target = 1.5,
      data_generator = mock_generator
    )
  )

  # Empty effect_name
  expect_error(
    brms_inla_power_sequential(
      formula = y ~ x, effect_name = character(0),
      effect_grid = c(0.5), sample_sizes = 20,
      metric = "direction", data_generator = mock_generator
    )
  )
})

test_that("Sequential engine runs and returns a classed result", {
  skip_if_not_installed("INLA")
  skip_on_cran()

  mock_generator <- function(n, effect) {
    effect_val <- if (is.list(effect)) effect[[1]] else effect
    data.frame(
      y = rnorm(n, effect_val, 1),
      treatment = rnorm(n, 0, 1)
    )
  }

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

  expect_s3_class(result, "brms_inla_power")
  expect_named(result, c("results", "summary", "settings"))
  expect_null(result$results)  # Sequential does not store per-sim results
  expect_s3_class(result$summary, "data.frame")
  expect_type(result$settings, "list")

  # Summary structure: conditional_power (renamed from 'assurance' in 1.3.0)
  expect_true("conditional_power" %in% names(result$summary))
  expect_false("assurance" %in% names(result$summary))
  expect_true(all(c("n", "sims_used", "power_direction") %in%
                    names(result$summary)))
  expect_true(all(result$summary$conditional_power >= 0 &
                    result$summary$conditional_power <= 1, na.rm = TRUE))
  expect_true(all(result$summary$sims_used <= 15))

  # print method now applies
  expect_output(print(result), "powerbrmsINLA")
})

test_that("Sequential engine handles multi-effect scenarios", {
  skip_if_not_installed("INLA")
  skip_on_cran()

  mock_generator_multi <- function(n, effect) {
    treatment_val <- if (is.list(effect)) effect$treatment else effect
    age_val <- if (is.list(effect)) effect$age else 0
    data.frame(
      y = rnorm(n, treatment_val + age_val, 1),
      treatment = rnorm(n, 0, 1),
      age = rnorm(n, 0, 1)
    )
  }

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

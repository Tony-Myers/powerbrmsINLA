# tests/testthat/test-compute_assurance.R
# Tests for compute_assurance() and supporting helpers

library(testthat)

# ---- Shared synthetic power_result (no INLA required) --------------------

make_syn_result <- function() {
  syn_summary <- data.frame(
    n               = rep(c(50, 100, 200), each = 3),
    treatment       = rep(c(0.2, 0.5, 0.8), 3),
    power_direction = c(0.40, 0.65, 0.85,
                        0.60, 0.82, 0.95,
                        0.72, 0.90, 0.98),
    power_threshold = c(0.30, 0.55, 0.75,
                        0.50, 0.72, 0.88,
                        0.62, 0.80, 0.92),
    stringsAsFactors = FALSE
  )
  list(
    summary  = syn_summary,
    settings = list(effect_name = "treatment")
  )
}

# Effect values present in the grid
EFFECTS <- c(0.2, 0.5, 0.8)


# ==========================================================================
# (a) Uniform weights reproduce the simple mean of per-cell powers
# ==========================================================================

test_that("uniform weights reproduce simple per-cell mean", {
  res <- make_syn_result()
  w   <- c("0.2" = 1/3, "0.5" = 1/3, "0.8" = 1/3)
  out <- compute_assurance(res, prior_weights = w, metric = "direction")

  # For each n, expected assurance = mean(power_direction across the 3 effect values)
  s         <- res$summary
  for (ni in c(50, 100, 200)) {
    s_n        <- s[s$n == ni, ]
    expected   <- mean(s_n$power_direction)
    actual     <- out$assurance$assurance[out$assurance$sample_size == ni]
    expect_equal(actual, expected, tolerance = 1e-10,
                 label = paste("n =", ni))
  }
})


# ==========================================================================
# (b) Degenerate prior (all weight on one cell) reproduces that cell's power
# ==========================================================================

test_that("degenerate prior reproduces single-cell power", {
  res <- make_syn_result()
  s   <- res$summary

  # All weight on treatment = 0.5
  w <- c("0.2" = 0, "0.5" = 1, "0.8" = 0)
  out <- compute_assurance(res, prior_weights = w, metric = "direction")

  for (ni in c(50, 100, 200)) {
    expected <- s$power_direction[s$n == ni & s$treatment == 0.5]
    actual   <- out$assurance$assurance[out$assurance$sample_size == ni]
    expect_equal(actual, expected, tolerance = 1e-10,
                 label = paste("n =", ni))
  }
})

test_that("degenerate prior works for threshold metric", {
  res <- make_syn_result()
  s   <- res$summary

  w <- c("0.2" = 1, "0.5" = 0, "0.8" = 0)
  out <- compute_assurance(res, prior_weights = w, metric = "threshold")

  expected_50 <- s$power_threshold[s$n == 50 & s$treatment == 0.2]
  actual_50   <- out$assurance$assurance[out$assurance$sample_size == 50]
  expect_equal(actual_50, expected_50, tolerance = 1e-10)
})


# ==========================================================================
# (c) Input validation: weights not summing to 1 (within tolerance 0.01)
# ==========================================================================

test_that("weights not summing to 1 cause an error", {
  res <- make_syn_result()

  # Weights sum to 0.9 — difference from 1 is 0.1 > 0.01
  w_bad <- c("0.2" = 0.3, "0.5" = 0.3, "0.8" = 0.3)
  expect_error(
    compute_assurance(res, prior_weights = w_bad),
    regexp = "sum to 1"
  )
})

test_that("weights just outside tolerance (sum = 1.02) cause an error", {
  res <- make_syn_result()
  w_bad <- c("0.2" = 0.34, "0.5" = 0.34, "0.8" = 0.34)  # sum = 1.02
  expect_error(
    compute_assurance(res, prior_weights = w_bad),
    regexp = "sum to 1"
  )
})

test_that("weights close to 1 within tolerance pass", {
  res <- make_syn_result()
  # Intentionally off by 0.005 (within default tolerance 0.01)
  w <- c("0.2" = 0.332, "0.5" = 0.333, "0.8" = 0.335)  # sum = 1.000
  # Should not error
  expect_no_error(
    compute_assurance(res, prior_weights = w)
  )
})


# ==========================================================================
# (d) print method runs without error
# ==========================================================================

test_that("print method runs without error (named vector weights)", {
  res <- make_syn_result()
  w   <- c("0.2" = 1/3, "0.5" = 1/3, "0.8" = 1/3)
  out <- compute_assurance(res, prior_weights = w)
  expect_no_error(print(out))
  # Returns invisibly
  expect_identical(print(out), out)
})

test_that("print method runs without error (distribution prior)", {
  res <- make_syn_result()
  out <- compute_assurance(
    res,
    prior_weights = list(dist = "normal", mean = 0.5, sd = 0.2)
  )
  expect_no_error(print(out))
})


# ==========================================================================
# Additional: return-value structure
# ==========================================================================

test_that("compute_assurance returns a list with the expected structure", {
  res <- make_syn_result()
  w   <- c("0.2" = 1/3, "0.5" = 1/3, "0.8" = 1/3)
  out <- compute_assurance(res, prior_weights = w, metric = "direction")

  expect_s3_class(out, "powerbrmsINLA_assurance")
  expect_named(out, c("assurance", "metric", "power_col", "prior_spec",
                       "weights", "eff_cols"))
  expect_s3_class(out$assurance, "data.frame")
  expect_true("sample_size" %in% names(out$assurance))
  expect_true("assurance"   %in% names(out$assurance))
  expect_equal(nrow(out$assurance), 3L)       # 3 sample sizes
  expect_true(all(out$assurance$assurance >= 0 & out$assurance$assurance <= 1))
  expect_equal(out$metric,    "direction")
  expect_equal(out$power_col, "power_direction")
  expect_equal(out$eff_cols,  "treatment")
})

test_that("compute_assurance increases with n (monotonicity check)", {
  res <- make_syn_result()
  w   <- c("0.2" = 0.2, "0.5" = 0.6, "0.8" = 0.2)
  out <- compute_assurance(res, prior_weights = w, metric = "direction")
  a   <- out$assurance$assurance
  # Assurance should not decrease as n increases (given the synthetic data)
  expect_true(a[1] <= a[2])
  expect_true(a[2] <= a[3])
})


# ==========================================================================
# Distribution-based prior_weights
# ==========================================================================

test_that("normal prior returns a valid assurance object", {
  res <- make_syn_result()
  out <- compute_assurance(
    res,
    prior_weights = list(dist = "normal", mean = 0.5, sd = 0.2)
  )
  expect_s3_class(out, "powerbrmsINLA_assurance")
  expect_true(all(is.finite(out$assurance$assurance)))
  expect_true(all(out$assurance$assurance >= 0 & out$assurance$assurance <= 1))
})

test_that("uniform prior gives same result as equal named weights", {
  res    <- make_syn_result()
  w_eq   <- c("0.2" = 1/3, "0.5" = 1/3, "0.8" = 1/3)
  out_eq <- compute_assurance(res, prior_weights = w_eq)
  out_un <- compute_assurance(res,
                              prior_weights = list(dist = "uniform"))
  expect_equal(out_eq$assurance$assurance,
               out_un$assurance$assurance,
               tolerance = 1e-10)
})

test_that("beta prior (mode/n) returns valid assurance", {
  res <- make_syn_result()
  out <- compute_assurance(
    res,
    prior_weights = list(dist = "beta", mode = 0.5, n = 5)
  )
  expect_true(all(is.finite(out$assurance$assurance)))
})


# ==========================================================================
# assurance_prior_weights() helper
# ==========================================================================

test_that("assurance_prior_weights returns named weights summing to 1", {
  effects <- c(0.2, 0.5, 0.8)
  w <- assurance_prior_weights(effects, dist = "normal", mean = 0.5, sd = 0.2)
  expect_named(w, as.character(effects))
  expect_equal(sum(w), 1, tolerance = 1e-10)
  expect_true(all(w >= 0))
})

test_that("assurance_prior_weights output is directly usable as prior_weights", {
  res     <- make_syn_result()
  effects <- c(0.2, 0.5, 0.8)
  w       <- assurance_prior_weights(effects, dist = "normal", mean = 0.5, sd = 0.2)
  expect_no_error(
    compute_assurance(res, prior_weights = w)
  )
})

test_that("assurance_prior_weights works for uniform and beta distributions", {
  effects <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  w_u <- assurance_prior_weights(effects, dist = "uniform")
  expect_equal(sum(w_u), 1, tolerance = 1e-10)

  w_b <- assurance_prior_weights(effects, dist = "beta",
                                  shape1 = 2, shape2 = 3)
  expect_equal(sum(w_b), 1, tolerance = 1e-10)
})


# ==========================================================================
# Edge cases and error handling
# ==========================================================================

test_that("power_result as plain data.frame is accepted", {
  s <- data.frame(
    n               = c(50, 50, 100, 100),
    treatment       = c(0.2, 0.8, 0.2, 0.8),
    power_direction = c(0.4, 0.8, 0.6, 0.9),
    stringsAsFactors = FALSE
  )
  w   <- c("0.2" = 0.5, "0.8" = 0.5)
  out <- compute_assurance(s, prior_weights = w)
  expect_s3_class(out, "powerbrmsINLA_assurance")
  expect_equal(nrow(out$assurance), 2L)
})

test_that("missing power column raises informative error", {
  res <- make_syn_result()
  w   <- c("0.2" = 1/3, "0.5" = 1/3, "0.8" = 1/3)
  expect_error(
    compute_assurance(res, prior_weights = w, metric = "rope"),
    regexp = "power_rope"
  )
})

test_that("distribution prior on multi-effect raises informative error", {
  # Multi-effect summary
  multi_summary <- data.frame(
    n               = rep(c(50, 100), each = 4),
    treatment       = rep(c(0.2, 0.2, 0.8, 0.8), 2),
    age             = rep(c(0.1, 0.3, 0.1, 0.3), 2),
    power_direction = runif(8, 0.4, 0.9),
    stringsAsFactors = FALSE
  )
  multi_result <- list(
    summary  = multi_summary,
    settings = list(effect_name = c("treatment", "age"))
  )
  expect_error(
    compute_assurance(multi_result,
                      prior_weights = list(dist = "normal", mean = 0.5, sd = 0.2)),
    regexp = "single-effect"
  )
})

test_that("multi-effect positional weights work", {
  multi_summary <- data.frame(
    n               = rep(c(50, 100), each = 4),
    treatment       = rep(c(0.2, 0.2, 0.8, 0.8), 2),
    age             = rep(c(0.1, 0.3, 0.1, 0.3), 2),
    power_direction = c(0.40, 0.45, 0.70, 0.75,
                        0.60, 0.65, 0.85, 0.90),
    stringsAsFactors = FALSE
  )
  multi_result <- list(
    summary  = multi_summary,
    settings = list(effect_name = c("treatment", "age"))
  )
  # 4 unique combinations; uniform weights
  w <- rep(0.25, 4)
  out <- compute_assurance(multi_result, prior_weights = w)
  expect_s3_class(out, "powerbrmsINLA_assurance")
  expect_equal(nrow(out$assurance), 2L)
  expect_true(all(out$assurance$assurance >= 0 & out$assurance$assurance <= 1))
})

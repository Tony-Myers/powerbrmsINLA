# tests/testthat/test-decide-sample-size.R
# Tests for decide_sample_size() — both assurance and conditional modes.

library(testthat)

# ---- Shared synthetic fixture -----------------------------------------------

make_syn <- function() {
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

# Normal design prior weights over the three effect values
make_weights <- function() {
  assurance_prior_weights(c(0.2, 0.5, 0.8), dist = "normal",
                          mean = 0.5, sd = 0.2)
}


# =============================================================================
# (a) Assurance mode: one row per metric, correct structure
# =============================================================================

test_that("assurance mode returns one row per requested metric", {
  res <- make_syn()
  w   <- make_weights()

  out <- decide_sample_size(res,
                             direction     = 0.80,
                             threshold     = 0.75,
                             prior_weights = w,
                             target        = 0.80)

  expect_s3_class(out, "powerbrmsINLA_sample_size")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  expect_setequal(out$metric, c("direction", "threshold"))
  expect_named(out, c("metric", "target", "n_recommended",
                       "assurance_achieved", "prior_description"))
})

test_that("assurance mode: single metric returns single row", {
  res <- make_syn()
  w   <- make_weights()

  out <- decide_sample_size(res,
                             direction     = 0.80,
                             prior_weights = w,
                             target        = 0.80)

  expect_equal(nrow(out), 1L)
  expect_equal(out$metric, "direction")
})

test_that("assurance mode: n_recommended is a valid sample size from the grid", {
  res <- make_syn()
  w   <- make_weights()

  out <- decide_sample_size(res,
                             direction     = 0.70,
                             prior_weights = w,
                             target        = 0.70)

  expect_false(is.na(out$n_recommended))
  expect_true(out$n_recommended %in% c(50, 100, 200))
  expect_true(out$assurance_achieved >= 0.70)
})

test_that("assurance mode: distribution prior (list) works", {
  res <- make_syn()
  dist_prior <- list(dist = "normal", mean = 0.5, sd = 0.2)

  out <- decide_sample_size(res,
                             direction     = 0.70,
                             prior_weights = dist_prior,
                             target        = 0.70)

  expect_s3_class(out, "powerbrmsINLA_sample_size")
  expect_equal(nrow(out), 1L)
  expect_false(is.na(out$n_recommended))
})

test_that("assurance mode: prior_description is populated", {
  res <- make_syn()
  w   <- make_weights()

  out <- decide_sample_size(res,
                             direction     = 0.70,
                             prior_weights = w,
                             target        = 0.70)

  expect_true(nchar(out$prior_description) > 0L)
})


# =============================================================================
# (b) Conditional mode: one row per effect-size value
# =============================================================================

test_that("conditional mode returns one row per effect-size value", {
  res <- make_syn()
  out <- decide_sample_size(res, direction = 0.80)

  expect_s3_class(out, "powerbrmsINLA_sample_size")
  # Three unique treatment values → three rows
  expect_equal(nrow(out), 3L)
  expect_true("treatment" %in% names(out))
  expect_true("n_recommended" %in% names(out))
})

test_that("conditional mode includes conditional power column", {
  res <- make_syn()
  out <- decide_sample_size(res, direction = 0.80)

  expect_true("cond_power_direction" %in% names(out))
  # Power at recommended n should be >= target for non-NA rows
  ok_rows <- !is.na(out$n_recommended)
  expect_true(all(out$cond_power_direction[ok_rows] >= 0.80))
})

test_that("conditional mode with threshold returns cond_power_threshold", {
  res <- make_syn()
  out <- decide_sample_size(res, threshold = 0.70)

  expect_true("cond_power_threshold" %in% names(out))
})

test_that("conditional mode with two metrics returns both power columns", {
  res <- make_syn()
  out <- decide_sample_size(res, direction = 0.60, threshold = 0.50)

  expect_true("cond_power_direction" %in% names(out))
  expect_true("cond_power_threshold" %in% names(out))
})


# =============================================================================
# (c) Backward compatibility: no prior_weights still works
# =============================================================================

test_that("calling without prior_weights uses conditional mode", {
  res <- make_syn()
  out <- decide_sample_size(res, direction = 0.80)

  expect_s3_class(out, "powerbrmsINLA_sample_size")
  expect_equal(attr(out, "mode"), "conditional")
  expect_false("metric" %in% names(out))  # assurance-mode column absent
})

test_that("targets list interface works in conditional mode", {
  res <- make_syn()
  out_direct <- decide_sample_size(res, direction = 0.80)
  out_list   <- decide_sample_size(res, targets = list(direction = 0.80))

  expect_equal(out_direct$n_recommended, out_list$n_recommended)
})

test_that("plain data.frame input works in conditional mode", {
  s <- data.frame(
    n               = rep(c(50, 100), each = 2),
    treatment       = rep(c(0.3, 0.7), 2),
    power_direction = c(0.55, 0.75, 0.78, 0.92),
    stringsAsFactors = FALSE
  )
  out <- decide_sample_size(s, direction = 0.75)

  expect_s3_class(out, "powerbrmsINLA_sample_size")
  expect_equal(nrow(out), 2L)
})

test_that("plain data.frame input works in assurance mode", {
  s <- data.frame(
    n               = rep(c(50, 100, 200), each = 3),
    treatment       = rep(c(0.2, 0.5, 0.8), 3),
    power_direction = c(0.40, 0.65, 0.85,
                        0.60, 0.82, 0.95,
                        0.72, 0.90, 0.98),
    stringsAsFactors = FALSE
  )
  # wrap as list so compute_assurance has settings
  res <- list(summary = s, settings = list(effect_name = "treatment"))
  w   <- make_weights()
  out <- decide_sample_size(res, direction = 0.70, prior_weights = w,
                             target = 0.70)

  expect_s3_class(out, "powerbrmsINLA_sample_size")
})


# =============================================================================
# (d) NA with informative message when no n meets the target
# =============================================================================

test_that("assurance mode returns NA when target is unachievable, with message", {
  res <- make_syn()
  w   <- make_weights()

  # Target of 0.999 is unachievable with this small grid
  expect_message(
    out <- decide_sample_size(res,
                               direction     = 0.999,
                               prior_weights = w,
                               target        = 0.999),
    regexp = "no sample size"
  )
  expect_equal(nrow(out), 1L)
  expect_true(is.na(out$n_recommended))
  expect_true(is.na(out$assurance_achieved))
})

test_that("conditional mode returns NA when target is unachievable, with message", {
  res <- make_syn()

  expect_message(
    out <- decide_sample_size(res, direction = 0.999),
    regexp = "no sample size met all targets"
  )
  expect_true(any(is.na(out$n_recommended)))
})


# =============================================================================
# (e) Sampled-SD aggregation in conditional mode
# =============================================================================

test_that("conditional mode aggregates across sampled_error_sd rows", {
  # Simulate a summary that has a sampled_error_sd column (variance uncertainty)
  s <- data.frame(
    n               = rep(c(50, 100), each = 4),
    treatment       = rep(c(0.3, 0.7), times = 4),
    sampled_error_sd = rep(c(0.9, 1.1), each = 4),
    power_direction  = c(0.55, 0.75, 0.55, 0.75,
                         0.80, 0.92, 0.80, 0.92),
    stringsAsFactors = FALSE
  )
  # Without aggregation this would yield 4 rows (2 treatments x 2 SD draws)
  # With aggregation it should yield 2 rows (one per treatment value)
  out <- decide_sample_size(s, direction = 0.75)

  expect_equal(nrow(out), 2L)
  expect_false("sampled_error_sd" %in% names(out))
})


# =============================================================================
# (f) print methods run without error
# =============================================================================

test_that("print method works in assurance mode", {
  res <- make_syn()
  w   <- make_weights()
  out <- decide_sample_size(res, direction = 0.70, prior_weights = w,
                             target = 0.70)
  expect_no_error(print(out))
  expect_identical(print(out), out)  # returns invisibly
})

test_that("print method works in conditional mode", {
  res <- make_syn()
  out <- decide_sample_size(res, direction = 0.80)
  expect_no_error(print(out))
  expect_identical(print(out), out)
})

test_that("print NA assurance row mentions criterion not achieved", {
  res <- make_syn()
  w   <- make_weights()
  suppressMessages(
    out <- decide_sample_size(res, direction = 0.999,
                               prior_weights = w, target = 0.999)
  )
  expect_output(print(out), regexp = "No sample size")
})

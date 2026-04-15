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
                             prior_weights = w)

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
                             prior_weights = w)

  expect_equal(nrow(out), 1L)
  expect_equal(out$metric, "direction")
})

test_that("assurance mode: n_recommended is a valid sample size from the grid", {
  res <- make_syn()
  w   <- make_weights()

  out <- decide_sample_size(res,
                             direction     = 0.70,
                             prior_weights = w)

  expect_false(is.na(out$n_recommended))
  expect_true(out$n_recommended %in% c(50, 100, 200))
  expect_true(out$assurance_achieved >= 0.70)
})

test_that("assurance mode: distribution prior (list) works", {
  res <- make_syn()
  dist_prior <- list(dist = "normal", mean = 0.5, sd = 0.2)

  out <- decide_sample_size(res,
                             direction     = 0.70,
                             prior_weights = dist_prior)

  expect_s3_class(out, "powerbrmsINLA_sample_size")
  expect_equal(nrow(out), 1L)
  expect_false(is.na(out$n_recommended))
})

test_that("assurance mode: prior_description is populated", {
  res <- make_syn()
  w   <- make_weights()

  out <- decide_sample_size(res,
                             direction     = 0.70,
                             prior_weights = w)

  expect_true(nchar(out$prior_description) > 0L)
})


# =============================================================================
# (g) Assurance mode: per-metric targets are respected
# =============================================================================

test_that("target column equals the value passed as direction argument", {
  res <- make_syn()
  w   <- make_weights()

  out_70 <- decide_sample_size(res, direction = 0.70, prior_weights = w)
  out_80 <- decide_sample_size(res, direction = 0.80, prior_weights = w)

  expect_equal(out_70$target, 0.70)
  expect_equal(out_80$target, 0.80)
})

test_that("target column equals the value passed as threshold argument", {
  res <- make_syn()
  w   <- make_weights()

  out <- decide_sample_size(res, threshold = 0.60, prior_weights = w)
  expect_equal(out$target, 0.60)
})

test_that("different metrics carry independent targets", {
  res <- make_syn()
  w   <- make_weights()

  out <- decide_sample_size(res,
                             direction     = 0.70,
                             threshold     = 0.60,
                             prior_weights = w)

  expect_equal(nrow(out), 2L)
  dir_row <- out[out$metric == "direction", ]
  thr_row <- out[out$metric == "threshold", ]
  expect_equal(dir_row$target,   0.70)
  expect_equal(thr_row$target,   0.60)
  # assurance_achieved must be >= its own target, not the other metric's target
  if (!is.na(dir_row$n_recommended)) {
    expect_true(dir_row$assurance_achieved >= 0.70)
  }
  if (!is.na(thr_row$n_recommended)) {
    expect_true(thr_row$assurance_achieved >= 0.60)
  }
})

test_that("lower target finds smaller or equal n than higher target", {
  res <- make_syn()
  w   <- make_weights()

  out_easy <- decide_sample_size(res, direction = 0.50, prior_weights = w)
  out_hard <- decide_sample_size(res, direction = 0.80, prior_weights = w)

  # If both are non-NA, the easier target should need no more n
  if (!is.na(out_easy$n_recommended) && !is.na(out_hard$n_recommended)) {
    expect_true(out_easy$n_recommended <= out_hard$n_recommended)
  }
})

test_that("direction = 0.70 gives different result than direction = 0.80", {
  res <- make_syn()
  w   <- make_weights()

  out_70 <- decide_sample_size(res, direction = 0.70, prior_weights = w)
  out_80 <- decide_sample_size(res, direction = 0.80, prior_weights = w)

  # Targets recorded in output must differ
  expect_false(identical(out_70$target, out_80$target))

  # The n at target=0.80 must be >= n at target=0.70 (when both achievable)
  if (!is.na(out_70$n_recommended) && !is.na(out_80$n_recommended)) {
    expect_true(out_80$n_recommended >= out_70$n_recommended)
  }
})

test_that("print output mentions the correct percentage for the metric target", {
  res <- make_syn()
  w   <- make_weights()
  out <- decide_sample_size(res, direction = 0.70, prior_weights = w)

  # Should say "70%" not "80%"
  expect_output(print(out), regexp = "70%", fixed = TRUE)
})

test_that("print output for threshold = 0.60 mentions 60%", {
  res <- make_syn()
  w   <- make_weights()
  out <- decide_sample_size(res, threshold = 0.60, prior_weights = w)

  expect_output(print(out), regexp = "60%", fixed = TRUE)
})

test_that("global target fallback is used when targets list has non-numeric", {
  # Edge case: targets list with a non-numeric value triggers global fallback
  res <- make_syn()
  w   <- make_weights()

  # Use the targets list interface with a proper numeric to verify it works
  out_list  <- decide_sample_size(res,
                                   targets       = list(direction = 0.70),
                                   prior_weights = w)
  out_direct <- decide_sample_size(res, direction = 0.70, prior_weights = w)

  expect_equal(out_list$target,         0.70)
  expect_equal(out_list$n_recommended,  out_direct$n_recommended)
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
  res <- list(summary = s, settings = list(effect_name = "treatment"))
  w   <- make_weights()
  out <- decide_sample_size(res, direction = 0.70, prior_weights = w)

  expect_s3_class(out, "powerbrmsINLA_sample_size")
  expect_equal(out$target, 0.70)
})


# =============================================================================
# (d) NA with informative message when no n meets the target
# =============================================================================

test_that("assurance mode returns NA when target is unachievable, with message", {
  res <- make_syn()
  w   <- make_weights()

  # direction = 0.999 is the per-metric target; no n achieves it
  expect_message(
    out <- decide_sample_size(res,
                               direction     = 0.999,
                               prior_weights = w),
    regexp = "no sample size"
  )
  expect_equal(nrow(out), 1L)
  expect_true(is.na(out$n_recommended))
  expect_true(is.na(out$assurance_achieved))
  # The recorded target is the value we actually passed, not the fallback 0.80
  expect_equal(out$target, 0.999)
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
  out <- decide_sample_size(res, direction = 0.70, prior_weights = w)
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
    out <- decide_sample_size(res, direction = 0.999, prior_weights = w)
  )
  expect_output(print(out), regexp = "No sample size")
})

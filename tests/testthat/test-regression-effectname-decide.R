# tests/testthat/test-regression-effectname-decide.R
#
# Regression tests for three fixes (all runnable without INLA):
#   #1  brms_inla_power() hard-errors on an effect_name that does not match a
#       formula-level fixed-effect term when the automatic data generator is in
#       use (and only warns when a custom data_generator is supplied).
#   #4  decide_sample_size() conditional mode requires at least one decision
#       target.
#   #5  decide_sample_size() does not mistake mean_sampled_* / sd_sampled_*
#       summary columns for effect-grid (design) columns.
#
# These tests deliberately avoid fitting any model: the validation in #1 runs
# before the INLA dependency check, and #4/#5 operate on plain data frames.

library(testthat)

# =============================================================================
# #1  effect_name validation
# =============================================================================

test_that("auto generator: non-matching effect_name is a hard error", {
  # Formula term is 'Condition'; the fitted coefficient level would be
  # 'Condition1'. Passing the coefficient level must fail, because the auto
  # generator keys on the formula term and would silently ignore it.
  expect_error(
    brms_inla_power(
      formula      = y ~ Condition + (1 | ID),
      effect_name  = "Condition1",
      effect_grid  = 1,
      sample_sizes = 40,
      nsims        = 5
    ),
    regexp = "must match a formula-level"
  )
})

test_that("auto generator: error fires before the INLA dependency check", {
  # The message must be the effect_name error, not the 'INLA is required'
  # error, regardless of whether INLA happens to be installed.
  err <- tryCatch(
    brms_inla_power(
      formula      = y ~ Condition + (1 | ID),
      effect_name  = "Nonsense",
      effect_grid  = 1,
      sample_sizes = 40,
      nsims        = 5
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("must match a formula-level", err))
  expect_false(grepl("INLA", err))
})

test_that("custom data_generator downgrades the mismatch to a warning", {
  # With a user-supplied generator we cannot know their naming convention, so a
  # mismatch should warn rather than stop. Execution then continues and stops
  # later (e.g. at the INLA check or inside the dummy generator); we swallow any
  # such downstream error and assert only that the warning was raised.
  dummy_gen <- function(n, effect) {
    data.frame(y = rnorm(n), Condition = factor(sample(c(0, 1), n, TRUE)),
               ID = factor(seq_len(n)))
  }
  expect_warning(
    tryCatch(
      brms_inla_power(
        formula        = y ~ Condition + (1 | ID),
        effect_name    = "Condition1",
        data_generator = dummy_gen,
        effect_grid    = 1,
        sample_sizes   = 40,
        nsims          = 5
      ),
      error = function(e) invisible(NULL)
    ),
    regexp = "must match a formula-level"
  )
})

test_that("a matching effect_name passes validation (no effect_name error)", {
  # We do not require INLA: if it is absent the call stops at the INLA check;
  # if it is present a (small) simulation may run. Either way the effect_name
  # validation must NOT be the cause of any error.
  outcome <- tryCatch(
    brms_inla_power(
      formula      = y ~ Condition + (1 | ID),
      effect_name  = "Condition",
      effect_grid  = 1,
      sample_sizes = 40,
      nsims        = 5
    ),
    error   = function(e) conditionMessage(e),
    warning = function(w) conditionMessage(w)
  )
  if (is.character(outcome)) {
    expect_false(grepl("must match a formula-level", outcome))
  } else {
    succeed("validation passed and the engine proceeded")
  }
})

# -----------------------------------------------------------------------------
# Root-cause check on the generator itself (the reason #1 matters): the effect
# is applied only when effect_name matches the formula term.
# -----------------------------------------------------------------------------

test_that("auto generator applies the effect only for a matching name", {
  gen_fn <- get(".auto_data_generator", asNamespace("powerbrmsINLA"))

  set.seed(1)
  gen_match <- gen_fn(formula = y ~ x, effect_name = "x",
                      error_sd = 0.001, group_sd = 0)
  d_match <- gen_match(2000, c(x = 1))
  slope_match <- stats::coef(stats::lm(y ~ x, data = d_match))[["x"]]
  # Effect of 1 requested -> estimated slope should be close to 1.
  expect_gt(slope_match, 0.9)
  expect_lt(slope_match, 1.1)

  set.seed(1)
  gen_miss <- gen_fn(formula = y ~ x, effect_name = "not_x",
                     error_sd = 0.001, group_sd = 0)
  # Request a deliberately huge effect (100). Because the name does not match
  # the term, the requested value is silently ignored and replaced by a small
  # noise coefficient ~ N(0, 0.3); the recovered slope therefore stays tiny and
  # nowhere near 100, whatever the RNG draw. This is the silent-failure that
  # fix #1 turns into a hard error at the engine level.
  d_miss <- gen_miss(2000, c(x = 100))
  slope_miss <- stats::coef(stats::lm(y ~ x, data = d_miss))[["x"]]
  expect_lt(abs(slope_miss), 5)
})

# =============================================================================
# #4  conditional mode requires at least one decision target
# =============================================================================

test_that("conditional mode errors when no decision target is supplied", {
  s <- data.frame(
    n               = c(50, 100),
    treatment       = c(0.5, 0.5),
    power_direction = c(0.60, 0.90),
    stringsAsFactors = FALSE
  )
  expect_error(
    decide_sample_size(s),
    regexp = "No decision target supplied"
  )
})

test_that("assurance mode still errors when no metric is supplied", {
  s <- data.frame(
    n               = rep(c(50, 100, 200), each = 3),
    treatment       = rep(c(0.2, 0.5, 0.8), 3),
    power_direction = c(0.40, 0.65, 0.85, 0.60, 0.82, 0.95, 0.72, 0.90, 0.98),
    stringsAsFactors = FALSE
  )
  res <- list(summary = s, settings = list(effect_name = "treatment"))
  w   <- assurance_prior_weights(c(0.2, 0.5, 0.8), dist = "normal",
                                 mean = 0.5, sd = 0.2)
  expect_error(
    decide_sample_size(res, prior_weights = w),
    regexp = "No metric targets supplied"
  )
})

test_that("supplying a single target still works (guard not over-eager)", {
  s <- data.frame(
    n               = c(50, 100),
    treatment       = c(0.5, 0.5),
    power_direction = c(0.60, 0.90),
    stringsAsFactors = FALSE
  )
  out <- decide_sample_size(s, direction = 0.80)
  expect_s3_class(out, "powerbrmsINLA_sample_size")
})

# =============================================================================
# #5  mean_sampled_* / sd_sampled_* are not treated as effect-grid columns
# =============================================================================

test_that("mean_sampled_/sd_sampled_ columns are excluded from effect grid", {
  # Real engine summaries carry per-cell SD moments named mean_sampled_error_sd
  # and sd_sampled_error_sd. Only 'treatment' is a genuine design factor here,
  # so the result must have one row per treatment value (2), not one row per
  # (treatment x SD-moment) combination (4).
  s <- data.frame(
    n                     = rep(c(50, 100), each = 2),
    treatment             = rep(c(0.3, 0.7), times = 2),
    mean_sampled_error_sd = rep(c(0.95, 1.05), each = 2),
    sd_sampled_error_sd   = rep(c(0.10, 0.12), each = 2),
    power_direction       = c(0.55, 0.80, 0.75, 0.92),
    stringsAsFactors      = FALSE
  )
  out <- decide_sample_size(s, direction = 0.75)

  expect_equal(nrow(out), 2L)
  expect_true("treatment" %in% names(out))
  expect_false("mean_sampled_error_sd" %in% names(out))
  expect_false("sd_sampled_error_sd" %in% names(out))
})

test_that("group_sd SD-moment columns are likewise excluded", {
  s <- data.frame(
    n                     = rep(c(50, 100), each = 2),
    treatment             = rep(c(0.3, 0.7), times = 2),
    mean_sampled_group_sd = rep(c(0.40, 0.55), each = 2),
    sd_sampled_group_sd   = rep(c(0.05, 0.07), each = 2),
    power_direction       = c(0.55, 0.80, 0.75, 0.92),
    stringsAsFactors      = FALSE
  )
  out <- decide_sample_size(s, direction = 0.75)

  expect_equal(nrow(out), 2L)
  expect_false("mean_sampled_group_sd" %in% names(out))
  expect_false("sd_sampled_group_sd" %in% names(out))
})

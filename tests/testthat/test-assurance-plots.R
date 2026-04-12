# tests/testthat/test-assurance-plots.R
# Tests for plot_assurance_curve(), plot_power_assurance_overlay(),
# and plot_design_prior()

library(testthat)

# ---- Shared synthetic fixtures -----------------------------------------------

make_power_result <- function() {
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

make_assurance <- function(pr = make_power_result(), metric = "direction") {
  compute_assurance(
    pr,
    prior_weights = list(dist = "normal", mean = 0.5, sd = 0.15),
    metric        = metric
  )
}

EFFECT_GRID <- c(0.2, 0.5, 0.8)


# =============================================================================
# plot_assurance_curve()
# =============================================================================

test_that("plot_assurance_curve returns a ggplot for a single assurance object", {
  a <- make_assurance()
  p <- plot_assurance_curve(a)
  expect_s3_class(p, "ggplot")
})

test_that("plot_assurance_curve returns a ggplot for a list of assurance objects", {
  pr <- make_power_result()
  a1 <- compute_assurance(pr, list(dist = "normal", mean = 0.5, sd = 0.10))
  a2 <- compute_assurance(pr, list(dist = "normal", mean = 0.5, sd = 0.30))
  p  <- plot_assurance_curve(list(a1, a2), labels = c("Informative", "Diffuse"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_assurance_curve accepts target = NULL (no reference line)", {
  a <- make_assurance()
  p <- plot_assurance_curve(a, target = NULL)
  expect_s3_class(p, "ggplot")
})

test_that("plot_assurance_curve accepts custom target value", {
  a <- make_assurance()
  p <- plot_assurance_curve(a, target = 0.90)
  expect_s3_class(p, "ggplot")
})

test_that("plot_assurance_curve errors when labels length mismatches list length", {
  pr <- make_power_result()
  a1 <- make_assurance(pr)
  a2 <- make_assurance(pr)
  expect_error(
    plot_assurance_curve(list(a1, a2), labels = c("Only one label")),
    regexp = "same length"
  )
})

test_that("plot_assurance_curve errors on non-assurance input", {
  expect_error(
    plot_assurance_curve(list("not an assurance object")),
    regexp = "powerbrmsINLA_assurance"
  )
})

test_that("plot_assurance_curve accepts custom title and subtitle", {
  a <- make_assurance()
  p <- plot_assurance_curve(a, title = "My title", subtitle = "My subtitle")
  expect_s3_class(p, "ggplot")
})


# =============================================================================
# plot_power_assurance_overlay()
# =============================================================================

test_that("plot_power_assurance_overlay returns a ggplot", {
  pr <- make_power_result()
  a  <- make_assurance(pr)
  p  <- plot_power_assurance_overlay(pr, a)
  expect_s3_class(p, "ggplot")
})

test_that("plot_power_assurance_overlay works for 'threshold' metric", {
  pr <- make_power_result()
  a  <- compute_assurance(pr,
                           prior_weights = list(dist = "normal", mean = 0.5, sd = 0.15),
                           metric        = "threshold")
  p  <- plot_power_assurance_overlay(pr, a, metric = "threshold")
  expect_s3_class(p, "ggplot")
})

test_that("plot_power_assurance_overlay errors on missing power column", {
  pr <- make_power_result()
  a  <- make_assurance(pr)
  expect_error(
    plot_power_assurance_overlay(pr, a, metric = "rope"),
    regexp = "power_rope"
  )
})

test_that("plot_power_assurance_overlay errors when assurance_result has wrong class", {
  pr <- make_power_result()
  expect_error(
    plot_power_assurance_overlay(pr, "not an assurance object"),
    regexp = "powerbrmsINLA_assurance"
  )
})

test_that("plot_power_assurance_overlay accepts custom title and subtitle", {
  pr <- make_power_result()
  a  <- make_assurance(pr)
  p  <- plot_power_assurance_overlay(pr, a,
                                      title    = "Custom title",
                                      subtitle = "Custom subtitle")
  expect_s3_class(p, "ggplot")
})

test_that("plot_power_assurance_overlay works with multi-effect summary", {
  multi_summary <- data.frame(
    n               = rep(c(50, 100), each = 4),
    treatment       = rep(c(0.2, 0.2, 0.8, 0.8), 2),
    age             = rep(c(0.1, 0.3, 0.1, 0.3), 2),
    power_direction = c(0.40, 0.45, 0.70, 0.75,
                        0.60, 0.65, 0.85, 0.90),
    stringsAsFactors = FALSE
  )
  multi_pr <- list(
    summary  = multi_summary,
    settings = list(effect_name = c("treatment", "age"))
  )
  a <- compute_assurance(multi_pr, prior_weights = rep(0.25, 4))
  p <- plot_power_assurance_overlay(multi_pr, a)
  expect_s3_class(p, "ggplot")
})


# =============================================================================
# plot_design_prior()
# =============================================================================

test_that("plot_design_prior returns a ggplot for a normal prior", {
  p <- plot_design_prior(
    prior_spec  = list(dist = "normal", mean = 0.5, sd = 0.15),
    effect_grid = EFFECT_GRID
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_design_prior returns a ggplot for a uniform prior", {
  p <- plot_design_prior(
    prior_spec  = list(dist = "uniform", min = 0.1, max = 0.9),
    effect_grid = EFFECT_GRID
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_design_prior returns a ggplot for a beta prior (shape1/shape2)", {
  p <- plot_design_prior(
    prior_spec  = list(dist = "beta", shape1 = 2, shape2 = 2),
    effect_grid = EFFECT_GRID
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_design_prior returns a ggplot for a beta prior (mode/n)", {
  p <- plot_design_prior(
    prior_spec  = list(dist = "beta", mode = 0.5, n = 5),
    effect_grid = EFFECT_GRID
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_design_prior returns a ggplot for multiple priors", {
  p <- plot_design_prior(
    prior_spec  = list(
      list(dist = "normal", mean = 0.5, sd = 0.10),
      list(dist = "normal", mean = 0.5, sd = 0.30)
    ),
    effect_grid = EFFECT_GRID,
    labels      = c("Informative", "Diffuse")
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_design_prior errors when labels length mismatches spec length", {
  expect_error(
    plot_design_prior(
      prior_spec  = list(
        list(dist = "normal", mean = 0.5, sd = 0.10),
        list(dist = "normal", mean = 0.5, sd = 0.30)
      ),
      effect_grid = EFFECT_GRID,
      labels      = c("Only one")
    ),
    regexp = "same length"
  )
})

test_that("plot_design_prior errors for unknown distribution", {
  expect_error(
    plot_design_prior(
      prior_spec  = list(dist = "cauchy", location = 0, scale = 1),
      effect_grid = EFFECT_GRID
    ),
    regexp = "Unsupported distribution"
  )
})

test_that("plot_design_prior accepts custom title, subtitle, n_points", {
  p <- plot_design_prior(
    prior_spec  = list(dist = "normal", mean = 0.5, sd = 0.15),
    effect_grid = EFFECT_GRID,
    n_points    = 50L,
    title       = "My prior",
    subtitle    = "For illustration"
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_design_prior errors on non-numeric effect_grid", {
  expect_error(
    plot_design_prior(
      prior_spec  = list(dist = "normal", mean = 0.5, sd = 0.15),
      effect_grid = "abc"
    ),
    regexp = "numeric"
  )
})

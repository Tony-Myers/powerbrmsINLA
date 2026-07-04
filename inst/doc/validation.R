## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## ----tc1-classical------------------------------------------------------------
# Exact classical power via the non-central t distribution
tc1_classical <- power.t.test(
  n           = 64,
  delta       = 0.5,
  sd          = 1,
  sig.level   = 0.025,
  type        = "two.sample",
  alternative = "one.sided"
)$power

round(tc1_classical, 4)

## ----tc1-bayes, eval = FALSE--------------------------------------------------
#  library(powerbrmsINLA)
#  
#  # Custom data generator: binary treatment, Gaussian response
#  two_sample_gen <- function(sigma) {
#    function(n, effect) {
#      delta <- as.numeric(effect[[1L]])
#      treat <- rbinom(n, 1L, 0.5)
#      y     <- delta * treat + rnorm(n, 0, sigma)
#      data.frame(treatment = treat, y = y, stringsAsFactors = FALSE)
#    }
#  }
#  
#  set.seed(101)
#  res_tc1 <- brms_inla_power(
#    formula        = y ~ treatment,
#    effect_name    = "treatment",
#    effect_grid    = 0.5,
#    sample_sizes   = 128L,        # total N (≈ 64 per group)
#    nsims          = 2000L,
#    prob_threshold = 0.975,       # Bayesian equivalent of one-sided alpha = 0.025
#    error_sd       = 1.0,
#    data_generator = two_sample_gen(1.0),
#    seed           = 101L,
#    progress       = "none"
#  )
#  
#  cat(sprintf("Bayesian direction power: %.4f\n", res_tc1$summary$power_direction))
#  cat(sprintf("Classical power.t.test:   %.4f\n", tc1_classical))
#  cat(sprintf("Absolute difference:      %.4f\n",
#              abs(res_tc1$summary$power_direction - tc1_classical)))

## ----tc2-classical------------------------------------------------------------
chen_n       <- c(20L, 100L, 133L, 200L)
chen_classic <- vapply(chen_n, function(ng) {
  power.t.test(
    n           = ng,
    delta       = 2.26,
    sd          = 6.536,
    sig.level   = 0.025,
    type        = "two.sample",
    alternative = "one.sided"
  )$power
}, numeric(1L))

data.frame(n_per_group = chen_n, classical_power = round(chen_classic, 4))

## ----tc2-bayes, eval = FALSE--------------------------------------------------
#  n_total_chen <- 2L * chen_n   # c(40, 200, 266, 400) total observations
#  
#  res_tc2 <- brms_inla_power(
#    formula        = y ~ treatment,
#    effect_name    = "treatment",
#    effect_grid    = 2.26,
#    sample_sizes   = n_total_chen,
#    nsims          = 2000L,
#    prob_threshold = 0.975,
#    error_sd       = 6.536,
#    data_generator = two_sample_gen(6.536),
#    seed           = 102L,
#    progress       = "none"
#  )
#  
#  print(res_tc2)

## ----tc3-bayesassurance, eval = FALSE-----------------------------------------
#  if (requireNamespace("bayesassurance", quietly = TRUE)) {
#    ba_result <- bayesassurance::assurance_nd_na(
#      n     = c(30L, 60L, 100L),
#      n0    = 0.01,
#      delta = 0.5,
#      sigma = 1.0,
#      alpha = 0.025
#    )
#    print(ba_result)
#  }

## ----tc3-powerbrmsINLA, eval = FALSE------------------------------------------
#  effect_grid <- seq(-0.5, 1.5, by = 0.1)
#  
#  # One-sample generator: constant 'mu' predictor; its INLA coefficient = mean
#  one_sample_gen <- function(n, effect) {
#    mu <- as.numeric(effect[["mu"]])
#    data.frame(y = rnorm(n, mean = mu, sd = 1), mu = 1L,
#               stringsAsFactors = FALSE)
#  }
#  
#  res_tc3 <- brms_inla_power(
#    formula        = y ~ mu,
#    effect_name    = "mu",
#    effect_grid    = effect_grid,
#    sample_sizes   = c(30L, 60L, 100L),
#    nsims          = 500L,
#    prob_threshold = 0.975,
#    error_sd       = 1.0,
#    data_generator = one_sample_gen,
#    seed           = 201L,
#    progress       = "none"
#  )
#  
#  # Design prior weights: N(0.5, 0.5^2) evaluated over the effect grid
#  w <- assurance_prior_weights(effect_grid, dist = "normal", mean = 0.5, sd = 0.5)
#  assur_tc3 <- compute_assurance(res_tc3, prior_weights = w)
#  print(assur_tc3)


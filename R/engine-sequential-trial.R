# engine-sequential-trial.R
# Sequential *trial* simulation: interim analyses of accumulating data.
#
# NOTE: This is distinct from brms_inla_power_sequential(), which adaptively
# stops the Monte Carlo *simulation* for efficiency. This engine simulates
# trials in which the *data* accumulate and are analysed at interim looks,
# returning the operating characteristics of a sequential Bayesian design.
#
# DRAFT for review/testing - candidate for a future release.


#' Sample a true effect from a design-prior specification
#'
#' @param spec List with a `dist` element: `"normal"` (mean, sd),
#'   `"uniform"` (min, max) or `"beta"` (shape1, shape2).
#' @param n Number of draws.
#' @return Numeric vector of length `n`.
#' @keywords internal
.sample_design_prior <- function(spec, n = 1L) {
  dist <- spec[["dist"]]
  if (is.null(dist) || !is.character(dist) || length(dist) != 1L) {
    stop("Design-prior list must have a character `dist` element ",
         '("normal", "uniform", or "beta").', call. = FALSE)
  }
  if (identical(dist, "normal")) {
    mu    <- spec[["mean"]] %||% 0
    sigma <- spec[["sd"]]   %||% 1
    if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
      stop("Normal design prior requires sd > 0.", call. = FALSE)
    }
    stats::rnorm(n, mean = mu, sd = sigma)
  } else if (identical(dist, "uniform")) {
    lo <- spec[["min"]]; hi <- spec[["max"]]
    if (is.null(lo) || is.null(hi) || lo >= hi) {
      stop("Uniform design prior requires min < max.", call. = FALSE)
    }
    stats::runif(n, min = lo, max = hi)
  } else if (identical(dist, "beta")) {
    a <- spec[["shape1"]]; b <- spec[["shape2"]]
    if (is.null(a) || is.null(b) || a <= 0 || b <= 0) {
      stop("Beta design prior requires shape1 > 0 and shape2 > 0.", call. = FALSE)
    }
    stats::rbeta(n, shape1 = a, shape2 = b)
  } else {
    stop("Unsupported design-prior distribution '", dist,
         "'. Supported: 'normal', 'uniform', 'beta'.", call. = FALSE)
  }
}


#' Simulate a sequential Bayesian trial with interim stopping rules
#'
#' @description
#' Simulates trials in which data **accumulate** and are analysed with INLA at
#' a prespecified schedule of interim looks. At each look a posterior
#' probability criterion is evaluated and the trial stops for success,
#' futility, or (for the ROPE metric) practical equivalence; otherwise it
#' continues to the maximum sample size. The function returns the design's
#' operating characteristics: the probability of each decision, the
#' distribution of stopping times, the expected sample size, and the average
#' effect estimate at early success stops (a direct measure of optional
#' stopping exaggeration).
#'
#' This complements [brms_inla_power()] (fixed-n designs) and is **not** the
#' same as [brms_inla_power_sequential()], which adaptively stops the Monte
#' Carlo simulation itself rather than simulating sequential trials.
#'
#' Because Bayesian posterior quantities obey the likelihood principle, the
#' *interpretation* of the final posterior is unaffected by the stopping rule;
#' the purpose of this simulation is to quantify the *frequency* properties of
#' the design (false-go rate, expected sample size, estimate exaggeration),
#' which do depend on the stopping rule.
#'
#' @details
#' **Accrual model.** For each simulated trial one dataset of the maximum size
#' is generated, and the analysis at a look of size `n` uses its first `n`
#' rows, so the data genuinely accumulate across looks. With the automatic
#' data generator, rows are exchangeable given the random effects (grouping
#' factors cycle across rows), which corresponds to all clusters recruiting in
#' parallel while observations accumulate within them. For a different accrual
#' pattern (for example, whole clusters arriving one at a time), supply a
#' custom `data_generator` whose rows are ordered by accrual time.
#'
#' **Decision rules.** Let `pr` be the posterior probability of the
#' prespecified alternative at a look (computed, as elsewhere in the package,
#' from a Normal approximation to the INLA fixed-effect posterior):
#' \itemize{
#'   \item `metric = "direction"`: `pr = P(effect > 0 | data)` (or `< 0` when
#'     `alternative = "less"`).
#'   \item `metric = "threshold"`: `pr = P(effect > effect_threshold | data)`
#'     (or `<` when `alternative = "less"`).
#'   \item `metric = "rope"`: two probabilities are monitored; the trial stops
#'     for "success" when `P(outside ROPE) >= prob_success` and for
#'     "equivalence" when `P(inside ROPE) >= prob_success`.
#' }
#' For the direction and threshold metrics the trial stops for success when
#' `pr >= prob_success` and for futility when `pr <= prob_futility` (futility
#' stopping is disabled when `prob_futility = NULL`). `prob_success` and
#' `prob_futility` may be vectors of `length(looks)` to implement stricter
#' early thresholds (O'Brien-Fleming-like behaviour).
#'
#' **What to report.** `summary` contains, per scenario: `p_success`,
#' `p_futility`, `p_equivalence`, `p_inconclusive`, `expected_n`,
#' `p_stop_early`, and `mean_est_at_success` next to `mean_true_at_success`
#' (their difference estimates the exaggeration of effect estimates caused by
#' early stopping). Under `true_effect` containing 0, `p_success` is the
#' false-go rate of the design.
#'
#' @param design Optional `"powerbrmsINLA_seq_design"` object from
#'   [sequential_design()]. When supplied, it overrides `formula`, `family`,
#'   `priors`, `effect_name`, `metric`, `alternative`, `looks`,
#'   `prob_success`, `prob_futility`, `effect_threshold`, `rope_bounds`, and
#'   `credible_level`, so the simulated operating characteristics correspond
#'   exactly to the design that will be monitored with
#'   [sequential_analysis()].
#' @param formula Model formula (brms syntax; random effects supported).
#' @param family Response family, e.g. `gaussian()`, `binomial()`.
#' @param family_control Optional INLA `control.family` list; takes precedence
#'   over translated priors, as in [brms_inla_power()].
#' @param priors Optional brms prior specification, translated via the
#'   built-in prior bridge (see [brms_inla_power()]).
#' @param data_generator Optional `function(n, effect)` returning a data frame
#'   with `n` rows **ordered by accrual time** (see Details). Defaults to the
#'   package's automatic generator.
#' @param effect_name Single character: the primary (monitored) effect.
#' @param true_effect Either a numeric vector of fixed true-effect scenarios
#'   (conditional sequential power per value), or a design-prior list such as
#'   `list(dist = "normal", mean = 0.5, sd = 0.15)`, in which case a fresh
#'   true effect is drawn for every simulated trial and the resulting success
#'   proportion is the design's *sequential assurance*.
#' @param looks Increasing integer vector of cumulative sample sizes (total
#'   observations, the same units as `sample_sizes` in [brms_inla_power()])
#'   at which interim analyses occur. The last element is the maximum n.
#' @param nsims Number of simulated trials per scenario.
#' @param metric `"direction"`, `"threshold"`, or `"rope"`.
#' @param alternative `"greater"` or `"less"`; the prespecified direction for
#'   the direction/threshold metrics. Ignored for `"rope"`.
#' @param prob_success Scalar or `length(looks)` vector of stop-for-success
#'   posterior probability thresholds (default 0.95).
#' @param prob_futility Scalar or `length(looks)` vector of stop-for-futility
#'   thresholds, or `NULL` (default) to disable futility stopping. Ignored
#'   with a warning for `metric = "rope"`.
#' @param effect_threshold Threshold for `metric = "threshold"` (on the model
#'   coefficient scale; use a negative value with `alternative = "less"`).
#' @param rope_bounds Length-2 numeric vector for `metric = "rope"`.
#' @param credible_level Credible level used for the reported interval widths.
#' @param error_sd,group_sd Scalar or distributional list, as in
#'   [brms_inla_power()]. With a distributional specification a fresh value is
#'   drawn once per simulated *trial* (a trial has one true nuisance value;
#'   it does not change between looks).
#' @param obs_per_group,predictor_means,predictor_sds Passed to the automatic
#'   data generator (see [brms_inla_power()]).
#' @param seed Integer seed.
#' @param inla_num_threads INLA threads ("outer:inner"); auto-detected if NULL.
#' @param family_args Named list of family-specific arguments for the
#'   automatic generator.
#' @param progress `"auto"`, `"text"`, or `"none"`.
#'
#' @return A list of class `"powerbrmsINLA_seq_trial"` with components:
#' \itemize{
#'   \item `results`: one row per simulated trial (scenario, true effect,
#'     decision, stopping n and look, posterior estimate at stopping).
#'   \item `summary`: operating characteristics per scenario (see Details).
#'   \item `look_summary`: per-look stopping distribution per scenario.
#'   \item `diagnostics`: failed/warned INLA fits per scenario.
#'   \item `settings`: a record of all design parameters.
#' }
#'
#' @examples
#' \dontrun{
#' looks <- c(40, 80, 120, 160, 200)
#'
#' # Conditional sequential power at fixed true effects (0 = false-go rate)
#' seq_fixed <- brms_inla_sequential_trial(
#'   formula          = y ~ treatment,
#'   effect_name      = "treatment",
#'   true_effect      = c(0, 0.5),
#'   looks            = looks,
#'   nsims            = 200,
#'   metric           = "threshold",
#'   effect_threshold = 0.2,
#'   prob_success     = 0.95,
#'   prob_futility    = 0.05,
#'   seed             = 123
#' )
#' print(seq_fixed)
#'
#' # Sequential assurance: draw the true effect from a design prior
#' seq_assur <- brms_inla_sequential_trial(
#'   formula          = y ~ treatment,
#'   effect_name      = "treatment",
#'   true_effect      = list(dist = "normal", mean = 0.5, sd = 0.15),
#'   looks            = looks,
#'   nsims            = 200,
#'   metric           = "threshold",
#'   effect_threshold = 0.2,
#'   prob_futility    = 0.05,
#'   seed             = 123
#' )
#' print(seq_assur)
#' }
#'
#' @export
brms_inla_sequential_trial <- function(
    design = NULL,
    formula = NULL,
    family = gaussian(),
    family_control = NULL,
    priors = NULL,
    data_generator = NULL,
    effect_name = NULL,
    true_effect = 0.5,
    looks = NULL,
    nsims = 200,
    metric = c("direction", "threshold", "rope"),
    alternative = c("greater", "less"),
    prob_success = 0.95,
    prob_futility = NULL,
    effect_threshold = 0,
    rope_bounds = NULL,
    credible_level = 0.95,
    error_sd = 1,
    group_sd = 0.5,
    obs_per_group = 10,
    predictor_means = NULL,
    predictor_sds = NULL,
    seed = 123,
    inla_num_threads = NULL,
    family_args = list(),
    progress = c("auto", "text", "none")
) {
  # ===== DESIGN OBJECT (optional): overrides the matching arguments =========
  if (!is.null(design)) {
    if (!inherits(design, "powerbrmsINLA_seq_design")) {
      stop("`design` must be a 'powerbrmsINLA_seq_design' object from ",
           "sequential_design().", call. = FALSE)
    }
    formula          <- design$formula
    family           <- design$family
    priors           <- design$priors
    effect_name      <- design$effect_name
    metric           <- design$metric
    alternative      <- design$alternative
    looks            <- design$looks
    prob_success     <- design$prob_success
    prob_futility    <- design$prob_futility
    effect_threshold <- design$effect_threshold
    rope_bounds      <- design$rope_bounds
    credible_level   <- design$credible_level
  }
  if (is.null(formula) || is.null(effect_name) || is.null(looks)) {
    stop("Supply either `design` (a sequential_design() object) or, at a ",
         "minimum, `formula`, `effect_name`, and `looks`.", call. = FALSE)
  }

  # ===== ARGUMENT VALIDATION (before dependency checks, so it is testable) ===
  metric      <- match.arg(metric)
  alternative <- match.arg(alternative)
  progress    <- match.arg(progress)

  stopifnot(is.character(effect_name))
  if (length(effect_name) != 1L) {
    stop("brms_inla_sequential_trial() monitors a single primary effect; ",
         "`effect_name` must have length 1. Additional fixed effects may ",
         "appear in the formula as nuisance terms.", call. = FALSE)
  }

  if (!is.numeric(looks) || length(looks) < 1L || any(!is.finite(looks)) ||
      any(looks <= 0) || is.unsorted(looks, strictly = TRUE)) {
    stop("`looks` must be a strictly increasing vector of positive sample ",
         "sizes (cumulative total observations at each interim analysis).",
         call. = FALSE)
  }
  looks <- as.integer(round(looks))
  n_looks <- length(looks)
  n_max   <- looks[n_looks]

  if (!is.numeric(prob_success) ||
      !(length(prob_success) %in% c(1L, n_looks)) ||
      any(prob_success <= 0.5) || any(prob_success >= 1)) {
    stop("`prob_success` must be a scalar or length(looks) vector with ",
         "values in (0.5, 1).", call. = FALSE)
  }
  prob_success <- rep_len(prob_success, n_looks)

  if (!is.null(prob_futility)) {
    if (identical(metric, "rope")) {
      warning("`prob_futility` is ignored for metric = 'rope'; equivalence ",
              "stopping uses `prob_success` applied to P(inside ROPE).",
              call. = FALSE)
      prob_futility <- NULL
    } else {
      if (!is.numeric(prob_futility) ||
          !(length(prob_futility) %in% c(1L, n_looks)) ||
          any(prob_futility <= 0) || any(prob_futility >= 0.5)) {
        stop("`prob_futility` must be NULL, or a scalar / length(looks) ",
             "vector with values in (0, 0.5).", call. = FALSE)
      }
      prob_futility <- rep_len(prob_futility, n_looks)
    }
  }

  if (identical(metric, "rope")) {
    if (is.null(rope_bounds) || length(rope_bounds) != 2L ||
        !is.numeric(rope_bounds) || rope_bounds[1L] >= rope_bounds[2L]) {
      stop("metric = 'rope' requires `rope_bounds` as a length-2 increasing ",
           "numeric vector, e.g. c(-0.2, 0.2).", call. = FALSE)
    }
  }

  validate_sd_spec(error_sd, "error_sd")
  validate_sd_spec(group_sd, "group_sd")
  error_sd_is_dist <- is.list(error_sd)
  group_sd_is_dist <- is.list(group_sd)

  # Scenarios: fixed numeric values, or a single design-prior scenario
  use_design_prior <- is.list(true_effect)
  if (use_design_prior) {
    # Fail fast on a bad specification (draw is discarded)
    invisible(.sample_design_prior(true_effect, 1L))
    scenarios <- "design_prior"
  } else {
    if (!is.numeric(true_effect) || length(true_effect) < 1L ||
        any(!is.finite(true_effect))) {
      stop("`true_effect` must be a finite numeric vector or a design-prior ",
           "list such as list(dist = 'normal', mean = 0.5, sd = 0.15).",
           call. = FALSE)
    }
    scenarios <- as.character(true_effect)
  }

  stopifnot(is.numeric(nsims), length(nsims) == 1L, nsims >= 1)
  nsims <- as.integer(nsims)

  # ===== DEPENDENCY CHECK ====================================================
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required for brms_inla_sequential_trial(). ",
         "See https://www.r-inla.org for installation instructions.",
         call. = FALSE)
  }

  set.seed(seed)

  # ===== AUTO-DETECT INLA THREADS ===========================================
  if (is.null(inla_num_threads)) {
    n_cores <- parallel::detectCores()
    if (!is.numeric(n_cores) || length(n_cores) != 1L || !is.finite(n_cores)) n_cores <- 1L
    inla_num_threads <- if (n_cores >= 4) "4:1" else if (n_cores >= 2) "2:1" else "1:1"
  }

  # ===== FAMILY / FORMULA / PRIORS ==========================================
  fam_map  <- .to_inla_family(family)
  fam_inla <- fam_map$inla
  needs_Ntrials <- fam_inla %in% c("binomial", "betabinomial")
  needs_E       <- fam_inla %in% c("poisson")

  prior_map <- .map_brms_priors_to_inla(
    priors,
    family_control_supplied = !is.null(family_control),
    inla_family = fam_inla
  )
  tf_alt <- .brms_to_inla_formula2(formula, hyper_by_re = prior_map$hyper_by_re)
  inla_formula_alt <- tf_alt$inla_formula
  re_specs         <- tf_alt$re_specs
  prior_map        <- .mark_unmatched_re_priors(prior_map, tf_alt$re_hyper_groups)
  prior_map        <- .audit_re_correlation_terms(prior_map, re_specs)

  # ===== DATA GENERATOR =====================================================
  using_auto_generator <- is.null(data_generator)
  if (using_auto_generator) {
    data_generator <- .auto_data_generator(
      formula = formula,
      effect_name = effect_name,
      family = family,
      family_args = family_args,
      error_sd = if (error_sd_is_dist) 1.0 else error_sd,
      group_sd = if (group_sd_is_dist) 0.5 else group_sd,
      obs_per_group = obs_per_group,
      predictor_means = predictor_means,
      predictor_sds = predictor_sds
    )
  } else {
    stopifnot(is.function(data_generator))
  }

  # ===== HELPERS ============================================================
  find_effect_in_fit <- function(eff, fitnames) {
    if (eff %in% fitnames) return(eff)
    candidates <- grep(paste0("^", eff), fitnames, value = TRUE)
    if (length(candidates) >= 1) return(candidates[1])
    NA_character_
  }

  fit_one_look <- function(dat) {
    inla_args <- list(
      formula = inla_formula_alt,
      data = dat,
      family = fam_inla,
      control.fixed = prior_map$control_fixed %||% list(),
      control.predictor = list(link = 1),
      control.family = family_control %||% prior_map$control_family %||% list(),
      num.threads = inla_num_threads,
      verbose = FALSE
    )
    if (needs_Ntrials && !is.null(dat$.Ntrials)) inla_args$Ntrials <- dat$.Ntrials
    if (needs_E && !is.null(dat$.E))            inla_args$E <- dat$.E
    if (!is.null(dat$.scale))                   inla_args$scale <- dat$.scale

    warns <- character(0L)
    fit <- tryCatch({
      withCallingHandlers(
        suppressMessages(do.call(INLA::inla, inla_args)),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    }, error = function(e) e)
    list(fit = fit, warned = length(warns) > 0L)
  }

  posterior_pr <- function(mu, sdv) {
    # Returns list(pr = main monitored probability, pr_in = P(inside ROPE))
    if (identical(metric, "direction")) {
      pr <- if (identical(alternative, "greater")) {
        1 - stats::pnorm(0, mu, sdv)
      } else {
        stats::pnorm(0, mu, sdv)
      }
      list(pr = pr, pr_in = NA_real_)
    } else if (identical(metric, "threshold")) {
      pr <- if (identical(alternative, "greater")) {
        1 - stats::pnorm(effect_threshold, mu, sdv)
      } else {
        stats::pnorm(effect_threshold, mu, sdv)
      }
      list(pr = pr, pr_in = NA_real_)
    } else { # rope
      p_in <- stats::pnorm(rope_bounds[2L], mu, sdv) -
        stats::pnorm(rope_bounds[1L], mu, sdv)
      list(pr = 1 - p_in, pr_in = p_in)
    }
  }

  show_progress <- (progress == "text") ||
    (progress == "auto" && interactive())

  # ===== SIMULATION =========================================================
  all_rows <- vector("list", length(scenarios))
  total_trials <- length(scenarios) * nsims
  done_trials  <- 0L

  for (sc_i in seq_along(scenarios)) {
    scen <- scenarios[sc_i]
    trial_rows <- vector("list", nsims)

    for (s in seq_len(nsims)) {
      # --- Per-trial truth -------------------------------------------------
      true_val <- if (use_design_prior) {
        .sample_design_prior(true_effect, 1L)
      } else {
        as.numeric(true_effect[sc_i])
      }

      cur_error_sd <- if (error_sd_is_dist) .sample_sd_spec(error_sd) else as.numeric(error_sd)
      cur_group_sd <- if (group_sd_is_dist) .sample_sd_spec(group_sd) else as.numeric(group_sd)
      if (using_auto_generator) {
        .gen_env <- environment(data_generator)
        if (error_sd_is_dist) .gen_env$error_sd <- cur_error_sd
        if (group_sd_is_dist) .gen_env$group_sd <- cur_group_sd
      }

      # --- Generate the FULL trial once; looks reuse its first n rows ------
      effects_named_vec <- stats::setNames(true_val, effect_name)
      dat_full <- data_generator(n_max, effects_named_vec)

      if (length(re_specs) > 0L) {
        for (re in re_specs) {
          gid <- as.integer(as.factor(dat_full[[re$group]]))
          if (isTRUE(re$has_intercept) && is.null(dat_full[[re$id_intercept]]))
            dat_full[[re$id_intercept]] <- gid
          if (!is.null(re$slope) && is.null(dat_full[[re$id_slope]]))
            dat_full[[re$id_slope]] <- gid
        }
      }

      # --- Interim looks ----------------------------------------------------
      decision   <- "inconclusive"
      stop_n     <- n_max
      stop_look  <- n_looks
      pr_stop    <- NA_real_
      est_mu     <- NA_real_
      est_sd     <- NA_real_
      ciw        <- NA_real_
      n_failed   <- 0L
      any_warn   <- FALSE
      any_fit_ok <- FALSE

      for (k in seq_len(n_looks)) {
        dat_k <- dat_full[seq_len(looks[k]), , drop = FALSE]
        fl <- fit_one_look(dat_k)
        any_warn <- any_warn || fl$warned
        fit <- fl$fit

        fitnames <- if (!inherits(fit, "error") && !is.null(fit$summary.fixed))
          rownames(fit$summary.fixed) else character()
        target <- find_effect_in_fit(effect_name, fitnames)

        if (inherits(fit, "error") || is.null(fit$summary.fixed) || is.na(target)) {
          n_failed <- n_failed + 1L
          next  # skip this look; the trial continues
        }
        any_fit_ok <- TRUE

        mu  <- as.numeric(fit$summary.fixed[target, "mean"])
        sdv <- as.numeric(fit$summary.fixed[target, "sd"])
        pp  <- posterior_pr(mu, sdv)

        z   <- stats::qnorm(1 - (1 - credible_level) / 2)
        est_mu <- mu; est_sd <- sdv; ciw <- 2 * z * sdv
        pr_stop <- pp$pr

        if (identical(metric, "rope")) {
          if (pp$pr >= prob_success[k]) {
            decision <- "success"; stop_n <- looks[k]; stop_look <- k; break
          }
          if (pp$pr_in >= prob_success[k]) {
            decision <- "equivalence"; stop_n <- looks[k]; stop_look <- k; break
          }
        } else {
          if (pp$pr >= prob_success[k]) {
            decision <- "success"; stop_n <- looks[k]; stop_look <- k; break
          }
          if (!is.null(prob_futility) && pp$pr <= prob_futility[k]) {
            decision <- "futility"; stop_n <- looks[k]; stop_look <- k; break
          }
        }
      }

      if (!any_fit_ok) decision <- "failed"

      trial_rows[[s]] <- tibble::tibble(
        scenario     = scen,
        sim          = s,
        true_effect  = true_val,
        decision     = decision,
        stop_n       = if (identical(decision, "failed")) NA_integer_ else stop_n,
        stop_look    = if (identical(decision, "failed")) NA_integer_ else stop_look,
        pr_at_stop   = pr_stop,
        est_mean     = est_mu,
        est_sd       = est_sd,
        ci_width     = ciw,
        looks_failed = n_failed,
        had_warning  = any_warn,
        sampled_error_sd = if (error_sd_is_dist) cur_error_sd else NA_real_,
        sampled_group_sd = if (group_sd_is_dist) cur_group_sd else NA_real_
      )

      done_trials <- done_trials + 1L
      if (show_progress && (done_trials %% max(1L, total_trials %/% 20L) == 0L)) {
        cat(sprintf("\rSequential trials: %d / %d (%d%%)",
                    done_trials, total_trials,
                    round(100 * done_trials / total_trials)))
        utils::flush.console()
      }
    }
    all_rows[[sc_i]] <- dplyr::bind_rows(trial_rows)
  }
  if (show_progress) cat("\n")

  res <- dplyr::bind_rows(all_rows)
  res_ok <- dplyr::filter(res, .data$decision != "failed")

  # ===== SUMMARIES ==========================================================
  # NOTE: do not create a summary column named `true_effect` before the
  # per-success means below — in summarise() a new column of the same name
  # masks the data column (including via .data), silently breaking the
  # subsetting. The mean is computed under a temporary name and renamed last.
  summ <- res_ok %>%
    dplyr::group_by(.data$scenario) %>%
    dplyr::summarise(
      mean_true_all  = mean(.data$true_effect),
      p_success      = mean(.data$decision == "success"),
      p_futility     = mean(.data$decision == "futility"),
      p_equivalence  = mean(.data$decision == "equivalence"),
      p_inconclusive = mean(.data$decision == "inconclusive"),
      expected_n     = mean(.data$stop_n),
      sd_n           = stats::sd(.data$stop_n),
      p_stop_early   = mean(.data$stop_n < n_max &
                              .data$decision %in% c("success", "futility", "equivalence")),
      mean_est_at_success  = if (any(.data$decision == "success"))
        mean(.data$est_mean[.data$decision == "success"], na.rm = TRUE) else NA_real_,
      mean_true_at_success = if (any(.data$decision == "success"))
        mean(.data$true_effect[.data$decision == "success"]) else NA_real_,
      nsims_ok       = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      true_effect = if (use_design_prior) NA_real_ else .data$mean_true_all
    ) %>%
    dplyr::select(dplyr::all_of(c(
      "scenario", "true_effect", "p_success", "p_futility", "p_equivalence",
      "p_inconclusive", "expected_n", "sd_n", "p_stop_early",
      "mean_est_at_success", "mean_true_at_success", "nsims_ok"
    )))

  scen_counts <- res_ok %>%
    dplyr::count(.data$scenario, name = "n_ok")

  look_summary <- res_ok %>%
    dplyr::filter(.data$decision %in% c("success", "futility", "equivalence")) %>%
    dplyr::count(.data$scenario, .data$stop_n, .data$decision) %>%
    dplyr::left_join(scen_counts, by = "scenario") %>%
    dplyr::mutate(prop = .data$n / .data$n_ok) %>%
    dplyr::arrange(.data$scenario, .data$stop_n) %>%
    dplyr::group_by(.data$scenario) %>%
    dplyr::mutate(cum_prop_stopped = cumsum(.data$prop)) %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::all_of(c("scenario", "stop_n", "decision",
                                  "n", "prop", "cum_prop_stopped")))

  diagnostics <- res %>%
    dplyr::group_by(.data$scenario) %>%
    dplyr::summarise(
      prop_warned        = mean(.data$had_warning),
      total_failed_looks = sum(.data$looks_failed),
      n_failed_trials    = sum(.data$decision == "failed"),
      .groups = "drop"
    )

  out <- list(
    results      = res,
    summary      = summ,
    look_summary = look_summary,
    diagnostics  = diagnostics,
    settings = list(
      formula          = formula,
      inla_formula     = inla_formula_alt,
      inla_family      = fam_inla,
      effect_name      = effect_name,
      true_effect      = true_effect,
      use_design_prior = use_design_prior,
      looks            = looks,
      n_max            = n_max,
      nsims            = nsims,
      metric           = metric,
      alternative      = alternative,
      prob_success     = prob_success,
      prob_futility    = prob_futility,
      effect_threshold = effect_threshold,
      rope_bounds      = rope_bounds,
      credible_level   = credible_level,
      error_sd         = error_sd,
      group_sd         = group_sd,
      obs_per_group    = obs_per_group,
      seed             = seed,
      inla_num_threads = inla_num_threads,
      prior_translation = prior_map$prior_audit,
      accrual_note = paste(
        "Each simulated trial generates one dataset at n_max;",
        "look k analyses its first looks[k] rows, so data accumulate.",
        "With the automatic generator, rows are exchangeable given the",
        "random effects (clusters recruit in parallel)."
      )
    )
  )
  class(out) <- c("powerbrmsINLA_seq_trial", "list")
  out
}


#' Print method for sequential trial simulation results
#'
#' @param x A `"powerbrmsINLA_seq_trial"` object from
#'   [brms_inla_sequential_trial()].
#' @param digits Number of digits for printed proportions.
#' @param ... Unused; present for S3 compatibility.
#' @return `x`, invisibly.
#' @export
print.powerbrmsINLA_seq_trial <- function(x, digits = 3L, ...) {
  cfg <- x$settings %||% list()

  cat("Sequential Bayesian trial simulation (powerbrmsINLA)\n")
  cat("=====================================================\n")
  cat("Effect      :", cfg$effect_name, "\n")
  cat("Looks (n)   :", paste(cfg$looks, collapse = ", "), "\n")
  cat("Metric      :", cfg$metric,
      if (!identical(cfg$metric, "rope")) paste0("(alternative: ", cfg$alternative, ")") else "",
      "\n")
  cat("Stop rules  : success >=", paste(unique(cfg$prob_success), collapse = "/"),
      if (!is.null(cfg$prob_futility))
        paste("; futility <=", paste(unique(cfg$prob_futility), collapse = "/"))
      else "; no futility stopping",
      "\n")
  if (isTRUE(cfg$use_design_prior)) {
    te <- cfg$true_effect
    cat("True effect : drawn from design prior (",
        te$dist, paste(unlist(te[setdiff(names(te), "dist")]), collapse = ", "),
        ") -> success rate is sequential ASSURANCE\n")
  } else {
    cat("True effect : fixed at", paste(cfg$true_effect, collapse = ", "),
        "-> success rates are CONDITIONAL sequential power\n")
  }
  cat("Simulations :", cfg$nsims, "per scenario\n")

  if (!is.null(x$diagnostics) && nrow(x$diagnostics) > 0L) {
    cat(sprintf(
      "INLA diagnostics: %.1f%% of trials had warnings; %d failed look(s); %d unusable trial(s).\n",
      mean(x$diagnostics$prop_warned, na.rm = TRUE) * 100,
      as.integer(sum(x$diagnostics$total_failed_looks, na.rm = TRUE)),
      as.integer(sum(x$diagnostics$n_failed_trials, na.rm = TRUE))
    ))
  }

  if (!is.null(x$summary) && nrow(x$summary) > 0L) {
    cat("\nOperating characteristics:\n")
    print(as.data.frame(lapply(x$summary, function(col)
      if (is.numeric(col)) round(col, digits) else col)), row.names = FALSE)

    n_max <- cfg$n_max %||% max(cfg$looks)
    sav <- 100 * (1 - x$summary$expected_n / n_max)
    cat(sprintf(
      "\nExpected sample size saving vs fixed n_max = %d: %s\n",
      as.integer(n_max),
      paste(sprintf("%s: %.0f%%", x$summary$scenario, sav), collapse = "; ")
    ))
    if (any(x$summary$p_success > 0, na.rm = TRUE)) {
      cat("Note: compare mean_est_at_success with mean_true_at_success to gauge\n")
      cat("the exaggeration of estimates caused by early stopping.\n")
    }
  }

  invisible(x)
}

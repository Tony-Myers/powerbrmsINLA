#' Core Bayesian Assurance / Power Simulation (Modern, Multi-Effect Ready)
#'
#' Provides Bayesian power analysis and assurance calculation using INLA
#' (Integrated Nested Laplace Approximation) for efficient computation.
#' Implements simulation-based power analysis for generalized linear mixed
#' models with automatic threading optimization.
#'
#' @param formula Model formula.
#' @param family brms GLM family (e.g., gaussian(), binomial()).
#' @param family_control Optional list for INLA's control.family.
#' @param Ntrials Optional vector for binomial trials.
#' @param E Optional vector for Poisson exposure.
#' @param scale Optional vector scale parameter for INLA families.
#' @param priors Optional brms::prior specification.
#' @param data_generator Optional function(n, effect) returning a dataset.
#' @param effect_name Character vector of fixed effect names.
#' @param effect_grid Vector/data.frame of effect values (supports multi-effect).
#'   For single effects, use a numeric vector. For multiple effects, use a
#'   data.frame with column names matching effect_name.
#' @param sample_sizes Vector of sample sizes.
#' @param nsims Number of simulations per cell.
#' @param power_threshold Decision probability threshold for summary.
#' @param precision_target Optional credible interval width target.
#' @param prob_threshold Posterior probability threshold for decision rules.
#' @param effect_threshold Effect-size threshold.
#' @param credible_level Credible interval level (default 0.95).
#' @param rope_bounds Optional Region of Practical Equivalence bounds (length 2 vector).
#' @param error_sd Gaussian residual standard deviation.
#' @param group_sd Random effects standard deviation.
#' @param obs_per_group Observations per group.
#' @param predictor_means Optional named list of predictor means.
#' @param predictor_sds Optional named list of predictor standard deviations.
#' @param seed Random seed.
#' @param inla_hyper Optional INLA-specific hyperparameters.
#' @param compute_bayes_factor Logical, compute Bayes Factor if TRUE.
#' @param bf_cutoff Numeric Bayes-factor threshold for declaring a "hit" (default 10).
#' @param bf_method Character. "sd" = Savage-Dickey at 0 (requires proper Normal prior on
#'   the tested coefficient); "marglik" = marginal-likelihood Bayes factor via INLA by
#'   comparing full vs reduced model (slower).
#' @param inla_num_threads Character string specifying INLA threading (e.g., "4:1"
#'   for 4 threads). If NULL (default), automatically detects optimal setting:
#'   "4:1" for 4+ cores, "2:1" for 2-3 cores, "1:1" otherwise.
#' @param progress One of "auto", "text", or "none" for progress display.
#' @param family_args List of arguments for family-specific data generators.
#' @return List with results, summary, and settings.
#' @export
brms_inla_power <- function(
    formula,
    family = gaussian(),
    family_control = NULL,
    Ntrials = NULL,
    E = NULL,
    scale = NULL,
    priors = NULL,
    data_generator = NULL,
    effect_name,
    effect_grid = 0.5,
    sample_sizes = c(50, 100, 200, 400),
    nsims = 200,
    power_threshold = 0.8,
    precision_target = NULL,
    prob_threshold = 0.95,
    effect_threshold = 0,
    credible_level = 0.95,
    rope_bounds = NULL,
    error_sd = 1,
    group_sd = 0.5,
    obs_per_group = 10,
    predictor_means = NULL,
    predictor_sds = NULL,
    seed = 123,
    inla_hyper = NULL,
    compute_bayes_factor = FALSE,
    bf_method = c("sd","marglik"),
    bf_cutoff = 10,
    inla_num_threads = NULL,
    progress = c("auto", "text", "none"),
    family_args = list()
) {
  # ===== INTERNAL UTILITY FUNCTIONS =====
  .to_inla_family_internal <- function(family) {
    family_name <- "unknown"
    if (is.character(family) && length(family) > 0) {
      family_name <- tolower(family[1])
    } else if (is.function(family)) {
      tryCatch({
        fam_obj <- family()
        if (is.list(fam_obj) && "family" %in% names(fam_obj)) {
          family_name <- tolower(fam_obj$family)
        }
      }, error = function(e) family_name <<- "unknown")
    } else if (is.list(family) && "family" %in% names(family)) {
      family_name <- tolower(family$family)
    } else {
      tryCatch({
        family_name <- tolower(as.character(family)[1])
      }, error = function(e) family_name <<- "unknown")
    }
    inla_family <- switch(family_name,
                          "gaussian" = "gaussian",
                          "binomial" = "binomial",
                          "bernoulli" = "binomial",
                          "poisson"  = "poisson",
                          "student"  = "T",
                          "negbinomial" = "nbinomial",
                          "negative binomial" = "nbinomial",
                          "beta"     = "beta",
                          "betabinomial" = "betabinomial",
                          "beta_binomial" = "betabinomial",
                          "lognormal" = "lognormal",
                          "weibull"   = "weibull",
                          "exponential" = "exponential",
                          "skew_normal" = "sn",
                          "skew-normal" = "sn",
                          "von_mises" = "vm",
                          "von-mises" = "vm",
                          "gev"       = "gev",
                          "gen_extreme_value" = "gev",
                          family_name
    )
    list(inla = inla_family, brms = family_name)
  }
  
  # Helper: pick the first log marginal likelihood safely (vector/matrix)
  .first_log_mlik <- function(fit) {
    if (inherits(fit, "error") || is.null(fit$mlik)) return(NA_real_)
    m <- fit$mlik
    if (is.matrix(m)) {
      val <- m[1, 1]
    } else {
      v <- as.numeric(m)
      val <- if (length(v) >= 1) v[1] else NA_real_
    }
    as.numeric(val)
  }
  
  # ===== PARAMETER VALIDATION =====
  stopifnot(is.character(effect_name), length(effect_name) > 0, nchar(effect_name[1]) > 0)
  set.seed(seed)
  progress   <- match.arg(progress)
  bf_method  <- match.arg(bf_method)
  stopifnot(is.numeric(bf_cutoff), length(bf_cutoff) == 1L, is.finite(bf_cutoff), bf_cutoff > 0)
  
  # ===== EFFECT GRID VALIDATION =====
  if (is.data.frame(effect_grid)) {
    grid_names <- colnames(effect_grid)
    if (!all(effect_name %in% grid_names)) {
      missing_effects <- effect_name[!effect_name %in% grid_names]
      stop("effect_grid column names must match effect_name.\n",
           "Missing columns: ", paste(missing_effects, collapse = ", "),
           "\nHave: ", paste(grid_names, collapse = ", "),
           "\nNeed: ", paste(effect_name, collapse = ", "))
    }
    effect_grid <- effect_grid[, effect_name, drop = FALSE]
    non_numeric <- !sapply(effect_grid, is.numeric)
    if (any(non_numeric)) {
      bad_cols <- names(effect_grid)[non_numeric]
      stop("All effect_grid columns must be numeric. Non-numeric: ",
           paste(bad_cols, collapse = ", "))
    }
    message("Multi-effect grid detected with ", nrow(effect_grid), " effect combinations")
  } else {
    if (length(effect_name) > 1) {
      stop("For multiple effect names, effect_grid must be a data.frame.")
    }
    if (!is.numeric(effect_grid)) stop("effect_grid must be numeric when given as a vector.")
    message("Single effect analysis for '", effect_name, "' with ", length(effect_grid), " values")
  }
  
  # ===== FORMULA VALIDATION =====
  formula_terms <- attr(terms(formula), "term.labels")
  formula_fixed <- formula_terms[!grepl("\\|", formula_terms)]
  for (eff in effect_name) {
    if (!any(grepl(paste0("\\b", eff, "\\b"),
                   c(formula_fixed, attr(terms(formula), "term.labels"))))) {
      warning("Effect name '", eff, "' not found in fixed-effects terms of the formula. ",
              "This may cause issues in data generation or model fitting.")
    }
  }
  
  # ===== AUTO-DETECT INLA THREADS =====
  if (is.null(inla_num_threads)) {
    n_cores <- parallel::detectCores()
    inla_num_threads <- if (n_cores >= 4) "4:1" else if (n_cores >= 2) "2:1" else "1:1"
  }
  
  # ===== PRIORS =====
  if (is.null(priors) || length(priors) == 0L) {
    message("No priors specified - using brms default weakly informative priors.")
    priors <- brms::prior("")
  }
  
  # ===== FAMILY SETUP =====
  fam_map    <- .to_inla_family_internal(family)
  fam_inla   <- fam_map$inla
  needs_Ntrials <- fam_inla %in% c("binomial", "betabinomial")
  needs_E       <- fam_inla %in% c("poisson")
  
  # ===== DATA GENERATOR =====
  if (is.null(data_generator)) {
    data_generator <- .auto_data_generator(
      formula = formula,
      effect_name = effect_name,
      family = family,
      family_args = family_args,
      error_sd = error_sd,
      group_sd = group_sd,
      obs_per_group = obs_per_group,
      predictor_means = predictor_means,
      predictor_sds = predictor_sds
    )
  } else {
    stopifnot(is.function(data_generator))
  }
  
  # ===== FORMULA + PRIORS MAP =====
  tf_alt <- .brms_to_inla_formula2(formula)
  inla_formula_alt <- tf_alt$inla_formula
  re_specs         <- tf_alt$re_specs
  prior_map        <- .map_brms_priors_to_inla(priors)
  
  # ===== GUARDS =====
  if (!is.null(rope_bounds) && length(rope_bounds) == 1) {
    warning("`rope` should be length-2 (e.g., c(0.1, 0.3)). ",
            "Using a single bound (", rope_bounds, "). ",
            "Did you write c(0.1:0.3) by mistake?", call. = FALSE)
  }
  
  # ===== HELPERS =====
  find_effect_in_fit <- function(eff, fitnames) {
    if (eff %in% fitnames) return(eff)
    candidates <- grep(paste0("^", eff), fitnames, value = TRUE)
    if (length(candidates) >= 1) return(candidates[1])
    NA_character_
  }
  get_prior_for_coef <- function(eff, prior_map_mean, prior_map_prec) {
    if (!is.null(prior_map_mean[[eff]])) {
      mean_val <- prior_map_mean[[eff]]
      sd_val <- if (!is.null(prior_map_prec[[eff]]) && prior_map_prec[[eff]] > 0)
        sqrt(1 / prior_map_prec[[eff]]) else NA_real_
      return(list(mean = mean_val, sd = sd_val))
    }
    eff_base <- sub("^(.*?)[0-9]+$", "\\1", eff)
    if (!is.null(prior_map_mean[[eff_base]])) {
      mean_val <- prior_map_mean[[eff_base]]
      sd_val <- if (!is.null(prior_map_prec[[eff_base]]) && prior_map_prec[[eff_base]] > 0)
        sqrt(1 / prior_map_prec[[eff_base]]) else NA_real_
      return(list(mean = mean_val, sd = sd_val))
    }
    list(mean = NA_real_, sd = NA_real_)
  }
  
  # ===== SIMULATION SETUP =====
  is_multi     <- is.data.frame(effect_grid)
  effect_rows  <- if (is_multi) seq_len(nrow(effect_grid)) else effect_grid
  total_steps  <- length(sample_sizes) * length(effect_rows) * nsims
  show_progress <- progress %in% c("auto", "text") && interactive()
  step <- 0L
  .simple_progress_bar <- function(step, total, width = 30) {
    done  <- round(width * step / total)
    bar   <- paste0(rep("=", done), collapse = "")
    space <- paste0(rep(" ", width - done), collapse = "")
    pct   <- round(100 * step / total)
    cat(sprintf("\r[%s%s] %3d%%", bar, space, pct))
    if (step == total) cat("\n")
    flush.console()
  }
  
  res_list <- vector("list", length(sample_sizes) * length(effect_rows))
  idx <- 0L
  
  for (n in sample_sizes) {
    for (eff_idx in effect_rows) {
      sim_rows <- vector("list", nsims)
      for (s in seq_len(nsims)) {
        # ===== GENERATE DATA =====
        if (is_multi) {
          eff_row <- effect_grid[eff_idx, , drop = FALSE]
          effects_named_vec <- setNames(as.numeric(eff_row), colnames(eff_row))
          dat <- data_generator(n, effects_named_vec)
        } else {
          effects_named_vec <- setNames(eff_idx, effect_name)
          dat <- data_generator(n, effects_named_vec)
        }
        
        # Add RE indices if needed
        if (length(re_specs) > 0L) {
          for (re in re_specs) {
            gid <- as.integer(as.factor(dat[[re$group]]))
            if (isTRUE(re$has_intercept) && is.null(dat[[re$id_intercept]]))
              dat[[re$id_intercept]] <- gid
            if (!is.null(re$slope) && is.null(dat[[re$id_slope]]))
              dat[[re$id_slope]] <- gid
          }
        }
        
        # ===== FIT INLA =====
        inla_args <- list(
          formula = inla_formula_alt,
          data = dat,
          family = fam_inla,
          control.fixed = prior_map$control_fixed %||% list(),
          control.predictor = list(link = 1),
          control.family = family_control %||% list(),
          num.threads = inla_num_threads,
          verbose = FALSE
        )
        if (needs_Ntrials && !is.null(dat$.Ntrials)) inla_args$Ntrials <- dat$.Ntrials
        if (needs_E && !is.null(dat$.E))            inla_args$E <- dat$.E
        if (!is.null(dat$.scale))                   inla_args$scale <- dat$.scale
        
        # Ensure INLA computes marginal likelihood for marglik
        if (isTRUE(compute_bayes_factor) && identical(bf_method, "marglik")) {
          inla_args$control.compute <- utils::modifyList(inla_args$control.compute %||% list(), list(mlik = TRUE))
        }
        
        fit <- tryCatch({
          suppressWarnings(suppressMessages(do.call(INLA::inla, inla_args)))
        }, error = function(e) e)
        
        # ===== EXTRACT =====
        fitnames <- if (!inherits(fit, "error") && !is.null(fit$summary.fixed))
          rownames(fit$summary.fixed) else character()
        target_coefs <- sapply(effect_name, find_effect_in_fit, fitnames = fitnames)
        
        if (inherits(fit, "error") || is.null(fit$summary.fixed) || any(is.na(target_coefs))) {
          sim_rows[[s]] <- tibble::tibble(
            sim = s, n = n, ok = FALSE,
            !!!setNames(as.list(rep(NA_real_, length(effect_name))), effect_name),
            post_prob_direction = NA_real_,
            post_prob_threshold = NA_real_,
            post_prob_rope = NA_real_,
            ci_width = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
            bf10 = NA_real_, log10_bf10 = NA_real_
          )
        } else {
          mean_b_vec <- sapply(target_coefs, function(nm) as.numeric(fit$summary.fixed[nm, "mean"]))
          sd_b_vec   <- sapply(target_coefs, function(nm) as.numeric(fit$summary.fixed[nm, "sd"]))
          
          prior_info <- get_prior_for_coef(target_coefs[1], prior_map$control_fixed$mean, prior_map$control_fixed$prec)
          prior_mean <- prior_info$mean
          prior_sd   <- prior_info$sd
          
          if (all(c("0.025quant", "0.975quant") %in% colnames(fit$summary.fixed))) {
            ci_lower <- as.numeric(fit$summary.fixed[target_coefs[1], "0.025quant"])
            ci_upper <- as.numeric(fit$summary.fixed[target_coefs[1], "0.975quant"])
          } else {
            ci_lower <- stats::qnorm((1 - credible_level) / 2, mean_b_vec[1], sd_b_vec[1])
            ci_upper <- stats::qnorm(1 - (1 - credible_level) / 2, mean_b_vec[1], sd_b_vec[1])
          }
          ci_width <- ci_upper - ci_lower
          
          eff_val <- if (is_multi) effects_named_vec[[effect_name[1L]]] else eff_idx
          dir_sign <- ifelse(eff_val >= 0, 1, -1)
          post_prob_direction <- if (dir_sign >= 0)
            1 - stats::pnorm(0, mean_b_vec[1], sd_b_vec[1]) else stats::pnorm(0, mean_b_vec[1], sd_b_vec[1])
          post_prob_threshold <- if (dir_sign >= 0)
            1 - stats::pnorm(effect_threshold, mean_b_vec[1], sd_b_vec[1]) else stats::pnorm(effect_threshold, mean_b_vec[1], sd_b_vec[1])
          post_prob_rope <- if (!is.null(rope_bounds) && length(rope_bounds) == 2L) {
            stats::pnorm(rope_bounds[2L], mean_b_vec[1], sd_b_vec[1]) -
              stats::pnorm(rope_bounds[1L], mean_b_vec[1], sd_b_vec[1])
          } else NA_real_
          
          # ---- Bayes factor (two methods) ----
          bf10 <- NA_real_; log10_bf10 <- NA_real_
          if (isTRUE(compute_bayes_factor)) {
            if (identical(bf_method, "sd")) {
              # Savage–Dickey: requires Normal prior on primary coef
              if (is.finite(prior_sd) && prior_sd > 0) {
                d_post0 <- stats::dnorm(0, mean_b_vec[1], sd_b_vec[1])
                d_pri0  <- stats::dnorm(0, mean = ifelse(is.finite(prior_mean), prior_mean, 0), sd = prior_sd)
                if (is.finite(d_post0) && d_post0 > 0 && is.finite(d_pri0)) {
                  bf10 <- d_pri0 / d_post0
                  log10_bf10 <- log10(bf10)
                }
              } else {
                warning(
                  "BF method 'sd' requires a proper Normal prior on '", effect_name[1],
                  "'. Bayes factors will be NA. Consider bf_method='marglik'.",
                  call. = FALSE
                )
              }
              
            } else { # bf_method == "marglik"
              # Ensure H1 has mlik; refit if needed
              fit_h1 <- fit
              if (is.null(fit_h1$mlik)) {
                inla_args_h1 <- inla_args
                inla_args_h1$control.compute <- utils::modifyList(inla_args_h1$control.compute %||% list(), list(mlik = TRUE))
                fit_h1 <- tryCatch(
                  suppressWarnings(suppressMessages(do.call(INLA::inla, inla_args_h1))),
                  error = function(e) e
                )
              }
              
              # Build H0 by dropping the primary tested effect
              drop_term <- paste(". ~ . -", effect_name[1])
              f0 <- tryCatch(stats::update(inla_formula_alt, drop_term), error = function(e) NULL)
              if (!is.null(f0)) {
                inla_args0 <- inla_args
                inla_args0$formula <- f0
                inla_args0$control.compute <- utils::modifyList(inla_args0$control.compute %||% list(), list(mlik = TRUE))
                
                fit0 <- tryCatch(
                  suppressWarnings(suppressMessages(do.call(INLA::inla, inla_args0))),
                  error = function(e) e
                )
                
                l1 <- .first_log_mlik(fit_h1)
                l0 <- .first_log_mlik(fit0)
                if (is.finite(l1) && is.finite(l0)) {
                  log10_bf10 <- (l1 - l0) / log(10)
                  bf10 <- 10^log10_bf10
                }
              }
            }
          }
          
          sim_rows[[s]] <- tibble::tibble(
            sim = s, n = n, ok = TRUE,
            !!!setNames(as.list(mean_b_vec), paste0("mean_", effect_name)),
            !!!setNames(as.list(sd_b_vec),   paste0("sd_",   effect_name)),
            !!!effects_named_vec,
            post_prob_direction = post_prob_direction,
            post_prob_threshold = post_prob_threshold,
            post_prob_rope = post_prob_rope,
            ci_width = ci_width, ci_lower = ci_lower, ci_upper = ci_upper,
            bf10 = bf10, log10_bf10 = log10_bf10
          )
        }
        
        step <- step + 1L
        if (show_progress) .simple_progress_bar(step, total_steps)
      } # sim loop
      
      idx <- idx + 1L
      res_list[[idx]] <- dplyr::bind_rows(sim_rows)
    }
  }
  
  # ===== COMBINE + SUMMARY =====
  res <- dplyr::bind_rows(res_list)
  group_vars <- if (is_multi) c("n", colnames(effect_grid)) else c("n", effect_name)
  summ <- res %>%
    dplyr::filter(ok) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      power_direction = mean(post_prob_direction >= prob_threshold, na.rm = TRUE),
      power_threshold = mean(post_prob_threshold >= prob_threshold, na.rm = TRUE),
      power_rope      = if (!is.null(rope_bounds)) mean(post_prob_rope <= (1 - prob_threshold), na.rm = TRUE) else NA_real_,
      avg_ci_width    = mean(ci_width, na.rm = TRUE),
      ci_coverage     = if (!is.null(precision_target)) mean(ci_width <= precision_target, na.rm = TRUE) else NA_real_,
      bf_hit_3        = mean(bf10 >= 3,  na.rm = TRUE),
      bf_hit_10       = mean(bf10 >= bf_cutoff, na.rm = TRUE),
      mean_log10_bf   = mean(log10_bf10, na.rm = TRUE),
      nsims_ok        = sum(ok, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    results = res,
    summary = summ,
    settings = list(
      formula = formula,
      inla_formula = inla_formula_alt,
      inla_family = fam_inla,
      effect_name = effect_name,
      effect_grid = effect_grid,
      sample_sizes = sample_sizes,
      nsims = nsims,
      prob_threshold = prob_threshold,
      effect_threshold = effect_threshold,
      credible_level = credible_level,
      rope_bounds = rope_bounds,
      power_threshold = power_threshold,
      precision_target = precision_target,
      compute_bayes_factor = compute_bayes_factor,
      bf_method = bf_method,
      bf_cutoff = bf_cutoff,
      prior_for_effect = list(
        mean = prior_map$control_fixed$mean,
        sd   = lapply(prior_map$control_fixed$prec, function(p)
          if (is.numeric(p) && p > 0) sqrt(1/p) else NA_real_)
      )
    )
  )
}
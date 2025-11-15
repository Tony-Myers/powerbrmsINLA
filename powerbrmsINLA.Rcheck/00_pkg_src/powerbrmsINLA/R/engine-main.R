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
#' @param inla_num_threads Character string specifying INLA threading (e.g., "4:1" 
#'   for 4 threads). If NULL (default), automatically detects optimal setting: 
#'   "4:1" for 4+ cores, "2:1" for 2-3 cores, "1:1" otherwise.
#' @param progress One of "auto", "text", or "none" for progress display.
#' @param family_args List of arguments for family-specific data generators.
#' @return List with results, summary, and settings.
#' 
#' @examples
#' \donttest{
#' # Basic usage with automatic INLA threading
#' results <- brms_inla_power(
#'   formula = outcome ~ treatment,
#'   effect_name = "treatment",
#'   effect_grid = c(0.2, 0.5, 0.8),
#'   sample_sizes = c(50, 100, 200),
#'   nsims = 3
#' )
#' print(results$summary)
#'
#' # Manual INLA threading control
#' results <- brms_inla_power(
#'   formula = outcome ~ treatment,
#'   effect_name = "treatment",
#'   effect_grid = c(0.2, 0.5, 0.8),
#'   sample_sizes = c(50, 100, 200),
#'   inla_num_threads = "8:1",  # Use 8 threads for faster computation
#'   nsims = 3
#' )
#'
#' # Multi-effect design with threading
#' effect_grid <- expand.grid(
#'   treatment = c(0, 0.3, 0.6),
#'   age_effect = c(0, 0.2)
#' )
#' results <- brms_inla_power(
#'   formula = outcome ~ treatment + age_effect,
#'   effect_name = c("treatment", "age_effect"),
#'   effect_grid = effect_grid,
#'   sample_sizes = c(100, 200, 400),
#'   nsims = 3
#' )
#' print(results$summary)
#' }
#' # Quick parameter check (runs instantly)
#' formals(brms_inla_power)
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
    inla_num_threads = NULL,
    progress = c("auto", "text", "none"),
    family_args = list()
) {
  
  # ===== INTERNAL UTILITY FUNCTIONS =====
  
  # Helper function: Map brms family to INLA family
  .to_inla_family_internal <- function(family) {
    # Extract family name safely
    family_name <- "unknown"
    
    # Handle character input
    if (is.character(family) && length(family) > 0) {
      family_name <- tolower(family[1])
    }
    # Handle function input (like gaussian(), binomial())
    else if (is.function(family)) {
      tryCatch({
        fam_obj <- family()
        if (is.list(fam_obj) && "family" %in% names(fam_obj)) {
          family_name <- tolower(fam_obj$family)
        }
      }, error = function(e) {
        family_name <<- "unknown"
      })
    }
    # Handle list/family object input
    else if (is.list(family) && "family" %in% names(family)) {
      family_name <- tolower(family$family)
    }
    # Last resort: try to convert to character
    else {
      tryCatch({
        family_name <- tolower(as.character(family)[1])
      }, error = function(e) {
        family_name <<- "unknown"
      })
    }
    
    # Map to INLA family name
    inla_family <- switch(family_name,
                          "gaussian" = "gaussian",
                          "binomial" = "binomial",
                          "bernoulli" = "binomial",
                          "poisson" = "poisson",
                          "student" = "T",
                          "negbinomial" = "nbinomial",
                          "negative binomial" = "nbinomial",
                          "beta" = "beta",
                          "betabinomial" = "betabinomial",
                          "beta_binomial" = "betabinomial",
                          "lognormal" = "lognormal",
                          "weibull" = "weibull",
                          "exponential" = "exponential",
                          "skew_normal" = "sn",
                          "skew-normal" = "sn",
                          "von_mises" = "vm",
                          "von-mises" = "vm",
                          "gev" = "gev",
                          "gen_extreme_value" = "gev",
                          family_name)  # Default: return as-is
    
    return(list(inla = inla_family, brms = family_name))
  }
  
  # ===== PARAMETER VALIDATION =====
  stopifnot(is.character(effect_name), length(effect_name) > 0, nchar(effect_name[1]) > 0)
  set.seed(seed)
  progress <- match.arg(progress)
  
  # ===== EFFECT GRID VALIDATION =====
  
  # Validate effect_grid structure and names
  if (is.data.frame(effect_grid)) {
    # Multi-effect case: validate column names match effect_name
    grid_names <- colnames(effect_grid)
    
    if (!all(effect_name %in% grid_names)) {
      missing_effects <- effect_name[!effect_name %in% grid_names]
      stop(
        "effect_grid column names must match effect_name.\n",
        "Missing columns in effect_grid: ", paste(missing_effects, collapse = ", "), "\n",
        "effect_grid columns: ", paste(grid_names, collapse = ", "), "\n",
        "effect_name: ", paste(effect_name, collapse = ", "), "\n\n",
        "Example fix:\n",
        "effect_grid <- data.frame(\n",
        "  ", paste(effect_name, " = c(...)", collapse = ",\n  "), "\n",
        ")"
      )
    }
    
    # Reorder columns to match effect_name order
    effect_grid <- effect_grid[, effect_name, drop = FALSE]
    
    # Validate all columns are numeric
    non_numeric <- !sapply(effect_grid, is.numeric)
    if (any(non_numeric)) {
      bad_cols <- names(effect_grid)[non_numeric]
      stop(
        "All effect_grid columns must be numeric.\n",
        "Non-numeric columns: ", paste(bad_cols, collapse = ", ")
      )
    }
    
    message("Multi-effect grid detected with ", nrow(effect_grid), " effect combinations")
    
  } else {
    # Single effect case: validate format
    if (length(effect_name) > 1) {
      stop(
        "For multiple effect names, effect_grid must be a data.frame.\n",
        "You specified ", length(effect_name), " effect names: ", 
        paste(effect_name, collapse = ", "), "\n",
        "But effect_grid is not a data.frame.\n\n",
        "Example fix:\n",
        "effect_grid <- expand.grid(\n",
        "  ", paste(effect_name, " = c(...)", collapse = ",\n  "), "\n",
        ")"
      )
    }
    
    if (!is.numeric(effect_grid)) {
      stop("effect_grid must be numeric when specified as a vector")
    }
    
    message("Single effect analysis for '", effect_name, "' with ", 
            length(effect_grid), " effect values")
  }
  
  # ===== FORMULA VALIDATION =====
  
  # Extract formula terms to validate effect names exist
  formula_terms <- attr(terms(formula), "term.labels")
  # Remove random effects for validation
  formula_fixed <- formula_terms[!grepl("\\|", formula_terms)]
  
  # Check if effect names appear in formula (allowing for interactions)
  for (eff in effect_name) {
    # Check direct match or as part of interaction
    if (!any(grepl(paste0("\\b", eff, "\\b"), c(formula_fixed, attr(terms(formula), "term.labels"))))) {
      warning(
        "Effect name '", eff, "' not found in formula terms.\n",
        "Formula terms: ", paste(formula_fixed, collapse = ", "), "\n",
        "This may cause issues in data generation or model fitting."
      )
    }
  }
  
  # ===== AUTO-DETECT OPTIMAL INLA THREADS =====
  if (is.null(inla_num_threads)) {
    n_cores <- parallel::detectCores()
    inla_num_threads <- if (n_cores >= 4) "4:1" else if (n_cores >= 2) "2:1" else "1:1"
  }
  
  # ===== USE BRMS DEFAULT PRIORS IF NOT PROVIDED =====
  if (is.null(priors) || length(priors) == 0L) {
    message("No priors specified - using brms default weakly informative priors.")
    priors <- brms::prior("")
  }
  
  # ===== SETUP FAMILY AND DATA REQUIREMENTS =====
  fam_map <- .to_inla_family_internal(family)
  fam_inla <- fam_map$inla
  needs_Ntrials <- fam_inla %in% c("binomial", "betabinomial")
  needs_E <- fam_inla %in% c("poisson")
  
  # ===== SETUP DATA GENERATOR =====
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
  
  # ===== FORMULA CONVERSION AND PRIOR MAPPING =====
  tf_alt <- .brms_to_inla_formula2(formula)
  inla_formula_alt <- tf_alt$inla_formula
  re_specs <- tf_alt$re_specs
  prior_map <- .map_brms_priors_to_inla(priors)
  
  # ===== HELPER FUNCTIONS =====
  
  # Helper: map user effect names to actual coefficient names in INLA fit
  find_effect_in_fit <- function(eff, fitnames) {
    if (eff %in% fitnames) return(eff)
    candidates <- grep(paste0("^", eff), fitnames, value = TRUE)
    if (length(candidates) >= 1) return(candidates[1])
    NA_character_
  }
  
  # Helper: get prior for coefficient
  get_prior_for_coef <- function(eff, prior_map_mean, prior_map_prec) {
    if (!is.null(prior_map_mean[[eff]])) {
      mean_val <- prior_map_mean[[eff]]
      sd_val <- if (!is.null(prior_map_prec[[eff]]) && prior_map_prec[[eff]] > 0)
        sqrt(1 / prior_map_prec[[eff]])
      else NA_real_
      return(list(mean = mean_val, sd = sd_val))
    }
    eff_base <- sub("^(.*?)[0-9]+$", "\\1", eff)
    if (!is.null(prior_map_mean[[eff_base]])) {
      mean_val <- prior_map_mean[[eff_base]]
      sd_val <- if (!is.null(prior_map_prec[[eff_base]]) && prior_map_prec[[eff_base]] > 0)
        sqrt(1 / prior_map_prec[[eff_base]])
      else NA_real_
      return(list(mean = mean_val, sd = sd_val))
    }
    list(mean = NA_real_, sd = NA_real_)
  }
  
  # ===== SIMULATION SETUP =====
  
  # Multi-effect grid detection (now properly validated)
  is_multi <- is.data.frame(effect_grid)
  effect_rows <- if (is_multi) seq_len(nrow(effect_grid)) else effect_grid
  total_steps <- length(sample_sizes) * length(effect_rows) * nsims
  show_progress <- progress %in% c("auto", "text") && interactive()
  step <- 0L
  
  .simple_progress_bar <- function(step, total, width = 30) {
    done <- round(width * step / total)
    bar <- paste0(rep("=", done), collapse = "")
    space <- paste0(rep(" ", width - done), collapse = "")
    pct <- round(100 * step / total)
    cat(sprintf("\r[%s%s] %3d%%", bar, space, pct))
    if (step == total) cat("\n")
    flush.console()
  }
  
  # ===== MAIN SIMULATION LOOP =====
  
  res_list <- vector("list", length(sample_sizes) * length(effect_rows))
  idx <- 0L
  
  for (n in sample_sizes) {
    for (eff_idx in effect_rows) {
      sim_rows <- vector("list", nsims)
      for (s in seq_len(nsims)) {
        # ===== GENERATE DATA WITH PROPER EFFECT NAMING =====
        if (is_multi) {
          # Get the specific row of effects
          eff_row <- effect_grid[eff_idx, , drop = FALSE]
          # Convert to named vector with correct names (already validated above)
          effects_named_vec <- setNames(as.numeric(eff_row), colnames(eff_row))
          dat <- data_generator(n, effects_named_vec)
        } else {
          # Single effect case
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
        
        # ===== FIT INLA MODEL =====
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
        if (needs_E && !is.null(dat$.E)) inla_args$E <- dat$.E
        if (!is.null(dat$.scale)) inla_args$scale <- dat$.scale
        
        fit <- tryCatch({
          # Suppress INLA warnings and messages during simulation
          suppressWarnings(suppressMessages({
            do.call(INLA::inla, inla_args)
          }))
        }, error = function(e) e)
        
        # ===== EXTRACT RESULTS =====
        fitnames <- if (!inherits(fit, "error") && !is.null(fit$summary.fixed))
          rownames(fit$summary.fixed) else character()
        target_coefs <- sapply(effect_name, find_effect_in_fit, fitnames = fitnames)
        
        if (inherits(fit, "error") || is.null(fit$summary.fixed) || any(is.na(target_coefs))) {
          # Failed simulation
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
          # Successful simulation - extract all metrics
          mean_b_vec <- sapply(target_coefs, function(nm) as.numeric(fit$summary.fixed[nm, "mean"]))
          sd_b_vec <- sapply(target_coefs, function(nm) as.numeric(fit$summary.fixed[nm, "sd"]))
          
          # Get prior for primary effect (for BF10)
          prior_info <- get_prior_for_coef(target_coefs[1], prior_map$control_fixed$mean, prior_map$control_fixed$prec)
          prior_mean <- prior_info$mean
          prior_sd <- prior_info$sd
          
          # Credible interval for first effect (main diagnostic)
          if (all(c("0.025quant", "0.975quant") %in% colnames(fit$summary.fixed))) {
            ci_lower <- as.numeric(fit$summary.fixed[target_coefs[1], "0.025quant"])
            ci_upper <- as.numeric(fit$summary.fixed[target_coefs[1], "0.975quant"])
          } else {
            ci_lower <- stats::qnorm((1 - credible_level) / 2, mean_b_vec[1], sd_b_vec[1])
            ci_upper <- stats::qnorm(1 - (1 - credible_level) / 2, mean_b_vec[1], sd_b_vec[1])
          }
          ci_width <- ci_upper - ci_lower
          
          # Posterior probability calculations
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
          
          # Bayes factor calculation
          bf10 <- NA_real_; log10_bf10 <- NA_real_
          if (isTRUE(compute_bayes_factor) && is.finite(prior_sd) && prior_sd > 0) {
            d_post0 <- stats::dnorm(0, mean_b_vec[1], sd_b_vec[1])
            d_pri0 <- stats::dnorm(0, mean = ifelse(is.finite(prior_mean), prior_mean, 0), sd = prior_sd)
            if (is.finite(d_post0) && is.finite(d_pri0) && d_post0 > 0) {
              bf10 <- d_pri0 / d_post0
              log10_bf10 <- log10(bf10)
            }
          }
          
          # Store results
          sim_rows[[s]] <- tibble::tibble(
            sim = s, n = n, ok = TRUE,
            !!!setNames(as.list(mean_b_vec), paste0("mean_", effect_name)),
            !!!setNames(as.list(sd_b_vec), paste0("sd_", effect_name)),
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
      } # end simulation loop
      
      idx <- idx + 1L
      res_list[[idx]] <- dplyr::bind_rows(sim_rows)
    }
  }
  
  # ===== COMBINE RESULTS AND COMPUTE SUMMARY =====
  
  res <- dplyr::bind_rows(res_list)
  
  # Summarize by effect grid columns (multi or single)
  group_vars <- if (is_multi) c("n", colnames(effect_grid)) else c("n", effect_name)
  summ <- res %>%
    dplyr::filter(ok) %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::summarise(
      power_direction = mean(post_prob_direction >= prob_threshold, na.rm = TRUE),
      power_threshold = mean(post_prob_threshold >= prob_threshold, na.rm = TRUE),
      power_rope     = if (!is.null(rope_bounds)) mean(post_prob_rope <= (1 - prob_threshold), na.rm = TRUE) else NA_real_,
      avg_ci_width   = mean(ci_width, na.rm = TRUE),
      ci_coverage    = if (!is.null(precision_target)) mean(ci_width <= precision_target, na.rm = TRUE) else NA_real_,
      bf_hit_3       = mean(bf10 >= 3, na.rm = TRUE),
      bf_hit_10      = mean(bf10 >= 10, na.rm = TRUE),
      mean_log10_bf  = mean(log10_bf10, na.rm = TRUE),
      nsims_ok       = sum(ok, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ===== RETURN RESULTS =====
  
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
      prior_for_effect = list(mean = prior_map$control_fixed$mean, sd = prior_map$control_fixed$prec)
    )
  )
}
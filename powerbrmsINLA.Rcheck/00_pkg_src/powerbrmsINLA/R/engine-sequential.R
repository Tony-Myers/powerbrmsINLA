#' Sequential Bayesian Assurance Simulation Engine (Modern, Multi-Effect Ready)
#'
#' Simulates assurance sequentially in batches, stopping early per cell based on Wilson confidence intervals.
#'
#' @param formula brms-style model formula.
#' @param family GLM family (e.g., gaussian(), binomial()).
#' @param family_control Optional list for INLA's control.family.
#' @param Ntrials Optional vector of binomial trial counts (for binomial families).
#' @param E Optional vector of exposures (for Poisson families).
#' @param scale Optional numeric vector for scale parameter in INLA.
#' @param priors brms prior specification object.
#' @param data_generator Optional function(n, effect) to simulate data.
#' @param effect_name Character vector of fixed effects to assess.
#' @param effect_grid Data frame or vector of effect values.
#' @param sample_sizes Vector of sample sizes.
#' @param metric Character; one of "direction", "threshold", "rope", or "bf" for Bayesian decision metric.
#' @param target Target assurance value for stopping.
#' @param prob_threshold Posterior probability threshold for decision metrics.
#' @param effect_threshold Effect-size threshold.
#' @param rope_bounds Numeric length-2 vector defining ROPE.
#' @param credible_level Credible interval level for Bayesian inference.
#' @param compute_bayes_factor Logical; TRUE if metric is "bf".
#' @param error_sd Residual standard deviation.
#' @param group_sd Standard deviation of random effects.
#' @param obs_per_group Number of observations per group.
#' @param predictor_means Optional named list of predictor means.
#' @param predictor_sds Optional named list of predictor standard deviations.
#' @param seed Random seed.
#' @param batch_size Number of simulations per sequential look.
#' @param min_sims Minimum simulations before early stopping.
#' @param max_sims Maximum simulations per cell.
#' @param ci_conf Confidence level for Wilson confidence intervals.
#' @param margin Margin around target for early stopping decision.
#' @param num_threads Number of INLA threads as character string (e.g., "1:1").
#' @param family_args List of family-specific args passed to data generator.
#' @param progress Logical; if TRUE, show progress messages.
#' @return List containing summary per cell and simulation settings.
#' @export

brms_inla_power_sequential <- function(
    formula,
    family = gaussian(),
    family_control = NULL,
    Ntrials = NULL,
    E = NULL,
    scale = NULL,
    priors = NULL,
    data_generator = NULL,
    effect_name,
    effect_grid,
    sample_sizes,
    metric = c("direction","threshold","rope","bf"),
    target = 0.8,
    prob_threshold = 0.95,
    effect_threshold = 0,
    rope_bounds = NULL,
    credible_level = 0.95,
    compute_bayes_factor = FALSE,
    error_sd = 1,
    group_sd = 0.5,
    obs_per_group = 10,
    predictor_means = NULL,
    predictor_sds = NULL,
    seed = 1,
    batch_size = 20,
    min_sims = 40,
    max_sims = 600,
    ci_conf = 0.95,
    margin = 0.02,
    num_threads = "1:1",
    family_args = list(),
    progress = TRUE
) {
  set.seed(seed)
  metric <- match.arg(metric)

  if (is.null(data_generator)) {
    data_generator <- .auto_data_generator(
      formula = formula, effect_name = effect_name,
      family = family, family_args = family_args,
      error_sd = error_sd, group_sd = group_sd, obs_per_group = obs_per_group,
      predictor_means = predictor_means, predictor_sds = predictor_sds
    )
  } else stopifnot(is.function(data_generator))

  tf <- .brms_to_inla_formula2(formula)
  f_inla <- tf$inla_formula
  re_specs <- tf$re_specs
  fam_inla <- .to_inla_family(family)$inla
  needs_N <- fam_inla %in% c("binomial","betabinomial")
  needs_E <- fam_inla %in% c("poisson")
  prior_map <- .map_brms_priors_to_inla(priors)

  # Multi-effect grid detection
  is_multi <- is.data.frame(effect_grid)
  effect_rows <- if (is_multi) seq_len(nrow(effect_grid)) else effect_grid
  total_cells <- length(sample_sizes) * length(effect_rows)

  prior_mean <- NA_real_; prior_sd <- NA_real_;
  # For multi-effect: get primary effect's prior for BF
  eff_main <- effect_name[1]
  if (isTRUE(compute_bayes_factor) && metric == "bf"
      && !is.null(prior_map$control_fixed$mean)
      && !is.null(prior_map$control_fixed$prec)
      && !is.null(prior_map$control_fixed$mean[[eff_main]])
      && !is.null(prior_map$control_fixed$prec[[eff_main]])) {
    prior_mean <- as.numeric(prior_map$control_fixed$mean[[eff_main]])
    prior_sd   <- sqrt(1 / as.numeric(prior_map$control_fixed$prec[[eff_main]]))
  }

  if (progress) message("Sequential assurance over ", total_cells, " cells (\u2026)")
  out <- vector("list", total_cells)
  idx <- 0L

  for (n in sample_sizes) {
    for (eff_idx in effect_rows) {
      hits <- 0L; trials <- 0L; stopped <- FALSE
      # The effect value(s)
      if (is_multi) {
        eff_vec <- as.list(effect_grid[eff_idx, , drop=FALSE])
        effects_named_vec <- setNames(unlist(eff_vec, use.names=TRUE), names(eff_vec))
        eff_val_main <- effects_named_vec[[eff_main]]
      } else {
        effects_named_vec <- setNames(eff_idx, effect_name)
        eff_val_main <- eff_idx
      }
      # For direction/threshold: sign for main
      dir_sign <- ifelse(eff_val_main >= 0, 1, -1)
      sims_used <- 0L

      while (!stopped && trials < max_sims) {
        b <- min(batch_size, max_sims - trials)
        for (s in seq_len(b)) {
          dat <- data_generator(n, effects_named_vec)
          if (length(re_specs) > 0L) {
            for (re in re_specs) {
              gid <- as.integer(as.factor(dat[[re$group]]))
              if (isTRUE(re$has_intercept) && is.null(dat[[re$id_intercept]]))
                dat[[re$id_intercept]] <- gid
              if (!is.null(re$slope) && is.null(dat[[re$id_slope]]))
                dat[[re$id_slope]] <- gid
            }
          }

          Ntrials_vec <- if (!is.null(dat$.Ntrials)) dat$.Ntrials else NULL
          E_vec      <- if (!is.null(dat$.E))      dat$.E else NULL
          scale_vec  <- if (!is.null(dat$.scale))  dat$.scale else NULL
          if (length(Ntrials_vec)==1L) Ntrials_vec <- rep(Ntrials_vec, n)
          if (length(E_vec)     ==1L) E_vec      <- rep(E_vec, n)
          if (length(scale_vec) ==1L) scale_vec  <- rep(scale_vec, n)

          fit <- tryCatch({
            inla_args <- list(
              formula          = f_inla,
              data             = dat,
              family           = fam_inla,
              control.fixed    = prior_map$control_fixed %||% list(),
              control.family   = family_control %||% list(),
              control.predictor= list(link=1),
              verbose          = FALSE,
              num.threads      = num_threads
            )
            if (needs_N && !is.null(Ntrials_vec)) inla_args$Ntrials <- Ntrials_vec
            if (needs_E && !is.null(E_vec)) inla_args$E <- E_vec
            if (!is.null(scale_vec)) inla_args$scale <- scale_vec
            do.call(INLA::inla, inla_args)
          }, error = function(e) e)
          # For multi-effect: try to retrieve each requested effect's coefficient
          fitnames <- if (!inherits(fit, "error") && !is.null(fit$summary.fixed))
            rownames(fit$summary.fixed) else character()
          target_coefs <- sapply(effect_name, function(eff) {
            if (eff %in% fitnames) return(eff)
            candidates <- grep(paste0("^", eff), fitnames, value = TRUE)
            if (length(candidates) >= 1) return(candidates[1])
            NA_character_
          })
          # Use first effect for metrics
          if (inherits(fit, "error") || is.null(fit$summary.fixed) || any(is.na(target_coefs))) next

          mean_b <- as.numeric(fit$summary.fixed[target_coefs[1], "mean"])
          sd_b   <- as.numeric(fit$summary.fixed[target_coefs[1], "sd"])
          success <- switch(metric,
                            "direction" = {
                              if (dir_sign >= 0) 1 - stats::pnorm(0, mean_b, sd_b) >= prob_threshold
                              else stats::pnorm(0, mean_b, sd_b) >= prob_threshold
                            },
                            "threshold" = {
                              thr <- effect_threshold
                              if (dir_sign >= 0) 1 - stats::pnorm(thr, mean_b, sd_b) >= prob_threshold
                              else stats::pnorm(thr, mean_b, sd_b) >= prob_threshold
                            },
                            "rope" = {
                              if (is.null(rope_bounds) || length(rope_bounds) != 2L) FALSE else {
                                p_in <- stats::pnorm(rope_bounds[2L], mean_b, sd_b) - stats::pnorm(rope_bounds[1L], mean_b, sd_b)
                                (1 - p_in) >= prob_threshold
                              }
                            },
                            "bf" = {
                              if (!is.finite(prior_sd) || prior_sd <= 0) FALSE else {
                                d_post0 <- stats::dnorm(0, mean_b, sd_b)
                                d_pri0  <- stats::dnorm(0, mean = ifelse(is.finite(prior_mean), prior_mean, 0), sd = prior_sd)
                                bf10 <- d_pri0 / d_post0
                                is.finite(bf10) && (bf10 >= prob_threshold)
                              }
                            }
          )

          hits <- hits + as.integer(success)
          trials <- trials + 1L
          sims_used <- trials
        }

        # Wilson CI-based early stop rule
        if (trials >= min_sims) {
          dec <- .should_stop_binom(hits, trials, target = target, margin = margin, conf = ci_conf)
          if (dec$stop) stopped <- TRUE
        }
      } # end while

      idx <- idx + 1L
      row_summary <- c(
        n = n,
        sims_used = sims_used,
        assurance = if (trials > 0) hits / trials else NA_real_,
        effect_val = eff_val_main,
        if (is_multi) effects_named_vec else NULL
      )
      # Set correct effect columns for summary output
      out[[idx]] <- as.data.frame(as.list(row_summary), stringsAsFactors=FALSE)
    }
  }

  # Combine results, set up column names
  res <- dplyr::bind_rows(out)
  # Construct summary table
  # multi-grid support: group/column names
  effect_cols <- if (is_multi) names(effect_grid) else "effect_val"
  summ <- res %>%
    dplyr::mutate(power_direction = if (metric == "direction") assurance else NA_real_,
                  power_threshold = if (metric == "threshold") assurance else NA_real_,
                  power_rope      = if (metric == "rope")      assurance else NA_real_,
                  bf_hit_10       = if (metric == "bf")         assurance else NA_real_) %>%
    dplyr::mutate(nsims_ok = sims_used) %>%
    dplyr::select(
      n, !!!effect_cols,
      power_direction, power_threshold, power_rope, bf_hit_10, nsims_ok,
      assurance, sims_used
    )

  list(
    results  = NULL,
    summary  = summ,
    settings = list(
      formula        = formula,
      inla_family    = fam_inla,
      effect_name    = effect_name,
      effect_grid    = effect_grid,
      sample_sizes   = sample_sizes,
      metric         = metric,
      target         = target,
      prob_threshold = prob_threshold,
      effect_threshold = effect_threshold,
      rope_bounds    = rope_bounds,
      credible_level = credible_level,
      compute_bayes_factor = compute_bayes_factor
    )
  )
}

#' Two-Stage Bayesian Assurance Simulation (Multi-Effect, User-Friendly API)
#'
#' Runs a two-stage Bayesian assurance simulation with formula-based multi-effect grids and adaptive refinement.
#'
#' @param formula Model formula.
#' @param effect_name Character vector of fixed effect names; must match formula terms.
#' @param effect_grid Data frame with columns named by effect_name specifying effect values.
#' @param n_range Numeric length-2 vector specifying sample size range.
#' @param stage1_k_n Number of grid points in stage 1.
#' @param stage1_nsims Number of simulations per cell in stage 1.
#' @param stage2_nsims Number of simulations per cell in stage 2.
#' @param refine_metric Metric used for refinement; one of "direction", "threshold", or "rope".
#' @param refine_target Target assurance for refined cells.
#' @param prob_threshold Posterior probability threshold for decision.
#' @param effect_threshold Effect-size threshold for decision metric.
#' @param obs_per_group Number of observations per group for grouping factors.
#' @param error_sd Residual standard deviation.
#' @param group_sd Standard deviation of random effects.
#' @param band Numeric width of the target refinement band.
#' @param expand Integer; how much to expand the refinement grid around candidates.
#' @param ... Additional arguments passed to internal functions.
#' @return A list with combined simulation results, summary, and stage parameters.
#' @export

brms_inla_power_two_stage <- function(
    formula,
    effect_name,            # Character vector, e.g. c("age", "treatment") - must match formula!
    effect_grid,            # data.frame with columns named by effect_name
    n_range,                # Numeric length-2: sample size range
    stage1_k_n = 8,
    stage1_nsims = 100,
    stage2_nsims = 400,
    refine_metric = c("direction", "threshold", "rope"),
    refine_target = 0.80,
    prob_threshold = 0.95,
    effect_threshold = 0.0,
    obs_per_group = NULL,
    error_sd = NULL,
    group_sd = 0.5,
    band = 0.06,
    expand = 1L,
    ...
) {
  refine_metric <- match.arg(refine_metric)
  stopifnot(is.data.frame(effect_grid), all(effect_name %in% colnames(effect_grid)))
  n_range <- sort(n_range)
  n1 <- unique(round(seq(n_range[1], n_range[2], length.out = stage1_k_n)))
  metric_col <- paste0("power_", refine_metric)

  # ----------- Stage 1: coarse grid over all effects -----------
  stage1 <- brms_inla_power(
    formula = formula,
    effect_name = effect_name,
    effect_grid = effect_grid,
    sample_sizes = n1,
    prob_threshold = prob_threshold,
    effect_threshold = effect_threshold,
    obs_per_group = obs_per_group,
    error_sd = error_sd,
    group_sd = group_sd,
    nsims = stage1_nsims,
    ...
  )

  # ----------- Find candidate cells near target band -----------
  lower <- refine_target - band
  upper <- refine_target + band
  summ1 <- stage1$summary
  cand <- dplyr::filter(summ1, .data[[metric_col]] >= lower, .data[[metric_col]] <= upper)
  if (nrow(cand) == 0L) {
    cand <- summ1 %>%
      dplyr::arrange(abs(.data[[metric_col]] - refine_target)) %>%
      dplyr::slice(seq_len(min(8L, nrow(.))))
}
  # ----------- Refine grid: expand around candidate cells ------
    expand_idx <- function(idx, maxlen) {
      as.integer(
        sort(
          unique(
            pmin(
              pmax(
                rep(idx, each = 2L * expand + 1L) + (-expand:expand),
                1L
              ),
              maxlen
            )
          )
        )
      )
    }

  effect_vals <- lapply(effect_name, function(col) {
    all_vals <- unique(effect_grid[[col]])
    # All candidate cell values as indices
    cand_idx <- match(cand[[col]], all_vals)
    expanded_idx <- sort(unique(unlist(lapply(cand_idx, expand_idx, maxlen = length(all_vals)))))
    all_vals[expanded_idx]
  })
  names(effect_vals) <- effect_name
  effect2_grid <- expand.grid(effect_vals)
  n2_idx <- match(cand$n, n1)
  n2 <- sort(unique(n1[unique(unlist(lapply(n2_idx, expand_idx, maxlen = length(n1))))]))

  # ----------- Stage 2: refined grid over expanded effects -----
  stage2 <- brms_inla_power(
    formula = formula,
    effect_name = effect_name,
    effect_grid = effect2_grid,
    sample_sizes = n2,
    prob_threshold = prob_threshold,
    effect_threshold = effect_threshold,
    obs_per_group = obs_per_group,
    error_sd = error_sd,
    group_sd = group_sd,
    nsims = stage2_nsims,
    ...
  )

  # ----------- Combine results and summary --------
  results_combined <- dplyr::bind_rows(stage1$results, stage2$results)
  summ <- dplyr::bind_rows(stage1$summary, stage2$summary) %>% dplyr::distinct()

  out <- list(
    results = results_combined,
    summary = summ,
    stage1 = list(effect_grid = effect_grid, grid_n = n1, nsims = stage1_nsims),
    stage2 = list(effect_grid = effect2_grid, grid_n = n2, nsims = stage2_nsims),
    params = list(
      refine_metric = refine_metric,
      refine_target = refine_target,
      prob_threshold = prob_threshold,
      effect_threshold = effect_threshold
    )
  )
  class(out) <- "brms_inla_power"
  out
}

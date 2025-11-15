#' Decide recommended sample size from power/assurance results
#'
#' Returns the smallest n per effect setting that meets user-specified targets.
#' Works with both brms_inla_power() and brms_inla_power_sequential() outputs.
#'
#' You can pass targets directly via arguments (direction, threshold, rope_in, bf10)
#' or via `targets = list(direction=..., threshold=..., rope_in=..., bf10=...)`.
#' Direct arguments take precedence if supplied.
#'
#' @param x A list with `$summary` (engine output) or a data.frame summary itself.
#' @param direction Numeric in \eqn{[0, 1]}, required power for `power_direction` (optional).
#' @param threshold Numeric in \eqn{[0, 1]}, required power for `power_threshold` (optional).
#' @param rope_in Numeric in \eqn{[0, 1]}, maximum allowed Pr(in ROPE) (optional).
#'   Note: since summaries usually contain `power_rope = Pr(outside ROPE) >= prob_threshold`,
#'   we compare `(1 - power_rope) <= rope_in` when `rope_in` is given.
#' @param bf10 Numeric Bayes-factor cutoff (e.g., 10). If provided, we look for a column
#'   named `bf_hit_<bf10>`; if not found, we fall back to any `bf_hit_*` column present.
#' @param bf_prop_min Numeric in \eqn{[0, 1]}, the minimum proportion of simulations that must
#'   achieve BF >= \code{bf10} (default 0).
#' @param targets Optional list alternative to the direct args. Ignored if any direct arg
#'   is non-NULL.
#' @return A data.frame with recommended n per effect combination and the rationale.
#' @export
decide_sample_size <- function(
    x,
    direction   = NULL,
    threshold   = NULL,
    rope_in     = NULL,
    bf10        = NULL,
    bf_prop_min = 0,
    targets     = NULL
) {
  # Get summary; accept list(summary = ..) or a data.frame
  s <- .get_summary_df(x)
  if (!"n" %in% names(s)) {
    stop("Summary data must contain a column 'n'.", call. = FALSE)
  }
  if (nrow(s) == 0L) {
    stop("Summary data has 0 rows; cannot decide sample size.", call. = FALSE)
  }
  s <- as.data.frame(s)
  s <- s[order(s$n), , drop = FALSE]
  
  # Merge explicit targets with list
  if (is.null(targets)) targets <- list()
  if (!is.null(direction)) targets$direction <- direction
  if (!is.null(threshold)) targets$threshold <- threshold
  if (!is.null(rope_in))   targets$rope_in   <- rope_in
  if (!is.null(bf10))      targets$bf10      <- bf10
  
  # Identify effect columns (everything except n and standard summary metrics)
  non_eff <- c(
    "n",
    "power_direction",
    "power_threshold",
    "power_rope",
    "avg_ci_width",
    "ci_coverage",
    "ci_width",
    "ci_lower",
    "ci_upper",
    "bf10_mean",
    "bf10_sd",
    "bf_hit_3",
    "bf_hit_10",
    grep("^bf_hit_", names(s), value = TRUE),
    "bf_median",
    "bf_min",
    "bf_max",
    "mean_log10_bf",
    "nsims_ok"
  )
  eff_cols <- setdiff(names(s), non_eff)
  
  ## --- Flexible BF proportion -------------------------------------------
  # We will create s$bf_prop if a BF target is requested:
  #  1. exact match to a bf_hit_<K> column -> use that
  #  2. otherwise, if x has $results with bf10 -> recompute Pr(BF10 >= bf10)
  #  3. otherwise, BF criterion is ignored with a message
  
  s$bf_prop <- NA_real_
  
  if (!is.null(targets$bf10)) {
    exact_col <- paste0("bf_hit_", as.integer(targets$bf10))
    
    if (exact_col %in% names(s)) {
      # Use precomputed proportion
      s$bf_prop <- s[[exact_col]]
    } else {
      # Try to recompute from per-simulation results
      r <- NULL
      if (is.list(x) && !is.null(x$results)) {
        r <- x$results
      }
      
      if (!is.null(r) && all(c("n", "bf10") %in% names(r))) {
        r <- as.data.frame(r)
        
        # Use the same effect columns that exist in both summary and results
        eff_cols_r <- intersect(eff_cols, names(r))
        group_cols_r <- c("n", eff_cols_r)
        
        bf_df <- r |>
          dplyr::group_by(dplyr::across(dplyr::all_of(group_cols_r))) |>
          dplyr::summarise(
            bf_prop = mean(bf10 >= targets$bf10, na.rm = TRUE),
            .groups = "drop"
          )
        
        join_cols <- c("n", eff_cols_r)
        s <- dplyr::left_join(
          s,
          bf_df,
          by = join_cols,
          suffix = c("", "_from_results")
        )
        
        # If bf_prop_from_results was created, move it into bf_prop
        if ("bf_prop_from_results" %in% names(s)) {
          s$bf_prop <- s$bf_prop_from_results
          s$bf_prop_from_results <- NULL
          message(
            "decide_sample_size(): recomputed Pr(BF10 >= ",
            targets$bf10,
            ") from per-simulation bf10 values."
          )
        } else {
          message(
            "decide_sample_size(): could not align per-simulation BF values with summary; ",
            "BF target will be ignored."
          )
        }
      } else {
        message(
          "decide_sample_size(): BF target requested but no matching bf_hit_* column ",
          "and no usable per-simulation bf10 in x$results; BF target will be ignored."
        )
      }
    }
  }
  
  ## --- Vectorised criteria -----------------------------------------------
  
  # direction
  ok_dir <- rep(TRUE, nrow(s))
  if (!is.null(targets$direction) && "power_direction" %in% names(s)) {
    ok_dir <- is.finite(s$power_direction) & (s$power_direction >= targets$direction)
  }
  
  # threshold
  ok_thr <- rep(TRUE, nrow(s))
  if (!is.null(targets$threshold) && "power_threshold" %in% names(s)) {
    ok_thr <- is.finite(s$power_threshold) & (s$power_threshold >= targets$threshold)
  }
  
  # ROPE: power_rope ~ Pr(outside ROPE); so Pr(in ROPE) = 1 - power_rope
  ok_rope <- rep(TRUE, nrow(s))
  if (!is.null(targets$rope_in) && "power_rope" %in% names(s)) {
    p_out <- s$power_rope
    p_in  <- ifelse(is.finite(p_out), 1 - p_out, NA_real_)
    ok_rope <- is.finite(p_in) & (p_in <= targets$rope_in)
  }
  
  # Bayes factor: use s$bf_prop if it is present and a BF target was requested
  ok_bf <- rep(TRUE, nrow(s))
  if (!is.null(targets$bf10) && "bf_prop" %in% names(s)) {
    ok_bf <- is.finite(s$bf_prop) & (s$bf_prop >= bf_prop_min)
  }
  
  ok_all <- ok_dir & ok_thr & ok_rope & ok_bf
  
  ## --- Group by effect grid and pick minimal n ---------------------------
  
  if (length(eff_cols) > 0L) {
    s$._group <- interaction(s[eff_cols], drop = TRUE, lex.order = TRUE)
    group_syms <- c("._group", eff_cols)
  } else {
    s$._group <- factor("all")
    group_syms <- "._group"
  }
  
  s$._ok_all <- ok_all
  
  out <- s |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_syms))) |>
    dplyr::summarise(
      n_recommended = if (any(._ok_all, na.rm = TRUE)) min(n[._ok_all], na.rm = TRUE) else NA_integer_,
      .groups       = "drop"
    )
  
  # Simple rationale
  out$rationale <- ifelse(
    is.na(out$n_recommended),
    "no n met all targets",
    "meets all"
  )
  
  # Put effect columns first in output
  if (length(eff_cols) > 0L) {
    out <- out[, c(eff_cols, "n_recommended", "rationale"), drop = FALSE]
  } else {
    out <- out[, c("n_recommended", "rationale"), drop = FALSE]
  }
  
  out
}

#' Plot decision/assurance curve across n
#'
#' @param x Engine result list (with $summary) or a data.frame.
#' @param y_metric One of "assurance","power_direction","power_threshold","power_rope","bf_hit_10".
#' @param target Optional horizontal target line.
#' @param effect_filter Named list for exact-match filtering (e.g., list(treatment=0.5)).
#' @param first_n_label If TRUE, annotate first n reaching target.
#' @return ggplot object.
#' @export
# Internal helper used for quick decision/assurance plots from a summary table
.plot_decision_assurance_curve_from_summary <- function(
    x,
    y_metric       = c("assurance","power_direction","power_threshold","power_rope","bf_hit_10"),
    target         = NULL,
    effect_filter  = NULL,
    first_n_label  = TRUE
) {
  y_metric <- match.arg(y_metric)
  df <- .get_summary_df(x)
  if (!"n" %in% names(df)) stop("Input must contain column 'n'.")
  
  # back-compat: "assurance" was often used as Pr(BF10 >= 10)
  if (identical(y_metric, "assurance")) y_metric <- "bf_hit_10"
  
  if (!is.null(effect_filter) && length(effect_filter)) {
    df <- .apply_effect_filters(df, effect_filter)
  }
  if (nrow(df) == 0L) stop("No rows left after filtering. Check effect_filter.")
  
  if (!y_metric %in% names(df)) stop("Column '", y_metric, "' not found in data.")
  
  df <- df[order(df$n), , drop = FALSE]
  
  n_first <- NA_real_
  if (!is.null(target) && is.finite(target)) {
    ok <- is.finite(df[[y_metric]]) & (df[[y_metric]] >= target)
    if (any(ok, na.rm = TRUE)) n_first <- min(df$n[ok], na.rm = TRUE)
  }
  
  can_draw_line <- isTRUE(dplyr::n_distinct(df$n) >= 2L)
  line_data <- if (can_draw_line) df else df[0, , drop = FALSE]
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n, y = .data[[y_metric]])) +
    ggplot2::geom_point() +
    ggplot2::geom_line(data = line_data, linewidth = 0.8) +
    ggplot2::labs(
      x     = "Sample size (n)",
      y     = y_metric,
      title = "Decision/assurance vs n"
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(target) && is.finite(target)) {
    p <- p + ggplot2::geom_hline(yintercept = target, linetype = "dashed")
  }
  if (first_n_label && is.finite(n_first)) {
    ymax <- max(df[[y_metric]], na.rm = TRUE)
    p <- p +
      ggplot2::geom_vline(xintercept = n_first, linetype = "dotted") +
      ggplot2::annotate(
        "text",
        x      = n_first,
        y      = ymax,
        label  = paste0("first n = ", n_first),
        vjust  = -0.5
      )
  }
  p
}
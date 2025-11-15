#' Bayes-factor assurance curve with Wilson CIs (multi-effect grid friendly)
#'
#' Plots the proportion of simulations in which BF10 meets or exceeds a threshold,
#' grouped by sample size and any effect grid variables.
#'
#' @param x Engine result (list with $results) or a data.frame with at least
#'   columns \code{n}, \code{bf10} and any effect columns (e.g., treatment).
#' @param cutoff Numeric Bayes factor threshold for a "hit" (default 10).
#' @param effect_filter Optional named list to filter effects,
#'   e.g. \code{list(treatment = 0.6)}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_bf_assurance_curve_smooth <- function(x, cutoff = 10, effect_filter = NULL) {
  df <- if (is.list(x) && !is.null(x$results)) x$results else x
  stopifnot(is.data.frame(df))
  if (!all(c("n", "bf10") %in% names(df))) {
    stop("Data must contain columns 'n' and 'bf10'.")
  }
  
  # effect columns = everything not obviously sim-level
  non_eff <- c(
    "sim", "n", "ok", "bf10", "log10_bf10",
    "ci_width", "ci_lower", "ci_upper",
    "post_prob_direction", "post_prob_threshold", "post_prob_rope",
    grep("^(mean_|sd_)", names(df), value = TRUE)
  )
  eff_cols <- setdiff(names(df), non_eff)
  
  # optional effect filter
  filtered_cols <- character(0)
  if (!is.null(effect_filter) && length(effect_filter)) {
    for (nm in names(effect_filter)) {
      if (!nm %in% names(df)) {
        stop("effect_filter refers to unknown column: ", nm)
      }
      df <- df[df[[nm]] == effect_filter[[nm]], , drop = FALSE]
      filtered_cols <- c(filtered_cols, nm)
    }
  }
  if (nrow(df) == 0L) stop("No rows left after filtering.")
  
  eff_cols <- setdiff(eff_cols, filtered_cols)
  
  # Wilson CI function (returns p_hat and CI for a binomial proportion)
  wilson <- function(h, t, conf = 0.95) {
    if (is.na(h) || is.na(t) || t == 0) {
      return(c(p = NA_real_, lo = NA_real_, hi = NA_real_))
    }
    z <- stats::qnorm(1 - (1 - conf) / 2)
    p <- h / t
    denom  <- 1 + z^2 / t
    centre <- (p + z^2 / (2 * t)) / denom
    half   <- (z * sqrt(p * (1 - p) / t + z^2 / (4 * t^2))) / denom
    c(p = p, lo = max(0, centre - half), hi = min(1, centre + half))
  }
  
  group_cols <- c("n", eff_cols)
  
  agg <- df %>%
    dplyr::mutate(hit = as.integer(is.finite(bf10) & bf10 >= cutoff)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      hits   = sum(hit, na.rm = TRUE),
      trials = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      p  = wilson(hits, trials)["p"],
      lo = wilson(hits, trials)["lo"],
      hi = wilson(hits, trials)["hi"]
    ) %>%
    dplyr::ungroup()
  
  # Build plot with proper grouping
  if (length(eff_cols) == 0L) {
    # No effect columns - simple plot
    p <- ggplot2::ggplot(agg, ggplot2::aes(x = n, y = p)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point()
  } else if (length(eff_cols) == 1L) {
    # One effect column - colour by it
    eff_var <- eff_cols[1L]
    p <- ggplot2::ggplot(
      agg,
      ggplot2::aes(
        x      = n,
        y      = p,
        group  = factor(!!rlang::sym(eff_var)),
        colour = factor(!!rlang::sym(eff_var))
      )
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = lo,
          ymax = hi,
          fill = factor(!!rlang::sym(eff_var))
        ),
        alpha  = 0.15,
        colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_viridis_d(name = eff_var) +
      ggplot2::scale_fill_viridis_d(name = eff_var)
  } else {
    # Multiple effect columns - colour by first, facet by rest
    eff_var    <- eff_cols[1L]
    facet_vars <- eff_cols[-1L]
    p <- ggplot2::ggplot(
      agg,
      ggplot2::aes(
        x      = n,
        y      = p,
        group  = factor(!!rlang::sym(eff_var)),
        colour = factor(!!rlang::sym(eff_var))
      )
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = lo,
          ymax = hi,
          fill = factor(!!rlang::sym(eff_var))
        ),
        alpha  = 0.15,
        colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_viridis_d(name = eff_var) +
      ggplot2::scale_fill_viridis_d(name = eff_var)
    
    if (length(facet_vars) == 1L) {
      p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_vars[1L])))
    } else {
      p <- p + ggplot2::facet_grid(
        as.formula(paste(paste(facet_vars, collapse = " + "), "~ ."))
      )
    }
  }
  
  p +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      title = paste0("Bayes-factor assurance (Pr(BF10 >= ", cutoff, "))"),
      x     = "Sample size (n)",
      y     = "Assurance"
    ) +
    ggplot2::theme_minimal()
}

#' Bayes-factor assurance curve (user-facing wrapper)
#'
#' This is the main function users should call to visualise
#' Bayes-factor "power" / assurance. It is a thin wrapper
#' around \code{plot_bf_assurance_curve_smooth()}, so all
#' existing behaviour is preserved.
#'
#' @param power_results Output from a brms_inla_power* function
#'   (or a data.frame with columns n, bf10 and any effect columns).
#' @param bf_threshold Numeric Bayes factor cutoff (default 10).
#' @param effect_filter Optional named list to filter effect-grid
#'   columns, e.g. \code{list(treatment = 0.3)}.
#'
#' @return A ggplot object.
#' @export
plot_bf_assurance_curve <- function(
    power_results,
    bf_threshold  = 10,
    effect_filter = NULL
) {
  plot_bf_assurance_curve_smooth(
    x             = power_results,
    cutoff        = bf_threshold,
    effect_filter = effect_filter
  )
}

#' Plot Expected Evidence (mean log10 BF10, Multi-Effect Grid Friendly)
#'
#' Plots the average log10 BF10 against any effect grid variable, grouped/faceted.
#'
#' @param power_results Simulation results from a `brms_inla_power*` function with `compute_bayes_factor = TRUE`.
#' @param x_effect Name of effect grid column for x-axis (default: first grid column).
#' @param facet_by Optional grid column(s) to facet by (default: NULL).
#' @param n Optional sample size to filter to (NULL means plot all; else one curve per grid/facet).
#' @param agg_fun Aggregation function if >1 entries per cell (default: mean).
#' @param title,subtitle Optional plot labels.
#' @return A ggplot object.
#' @export
plot_bf_expected_evidence <- function(
    power_results,
    x_effect  = NULL,
    facet_by  = NULL,
    n         = NULL,
    agg_fun   = mean,
    title     = NULL,
    subtitle  = NULL
) {
  s <- .get_summary_df(power_results)
  
  if (!"mean_log10_bf" %in% names(s)) {
    if ("bf10_log10_mean" %in% names(s)) {
      s$mean_log10_bf <- s$bf10_log10_mean
    } else {
      stop("Summary must contain 'mean_log10_bf' or 'bf10_log10_mean'.", call. = FALSE)
    }
  }
  
  # identify effect columns
  non_eff <- c("n", "mean_log10_bf", "bf10_log10_mean", "ci_width",
               "power_direction", "power_threshold", "power_rope")
  effect_names <- setdiff(names(s), non_eff)
  if (length(effect_names) == 0L) {
    stop("No effect grid columns found in summary.", call. = FALSE)
  }
  
  x_effect   <- x_effect %||% effect_names[1L]
  facet_names <- facet_by %||% character(0)
  facet_names <- intersect(facet_names, names(s))
  
  if (!x_effect %in% names(s)) {
    stop("x_effect column '", x_effect, "' not found in summary.", call. = FALSE)
  }
  
  # optional filter on n
  if (!is.null(n) && "n" %in% names(s)) {
    s <- s[s$n %in% n, , drop = FALSE]
  }
  
  group_cols <- c(x_effect, "n", facet_names)
  
  tmp <- s |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      mean_log10_bf = agg_fun(mean_log10_bf, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot2::ggplot(
    tmp,
    ggplot2::aes(x = !!rlang::sym(x_effect), y = mean_log10_bf)
  )
  
  if ("n" %in% names(tmp)) {
    p <- p +
      ggplot2::aes(
        group  = factor(n),
        colour = factor(n)
      ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_viridis_d(name = "n")
  } else {
    p <- p + ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point()
  }
  
  if (length(facet_names) == 1L) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_names[1L])))
  } else if (length(facet_names) > 1L) {
    p <- p + ggplot2::facet_grid(
      as.formula(paste(paste(facet_names, collapse = " + "), "~ ."))
    )
  }
  
  p +
    ggplot2::labs(
      x        = x_effect,
      y        = "Expected log10 BF10",
      title    = title %||% "Expected evidence (mean log10 BF10)",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
}
#' Plot Bayes Factor Heatmap (mean log10 BF10, Multi-Effect Grid Friendly)
#'
#' Heatmap of mean log10 BF10 as a function of two effect grid columns (x/y), with optional faceting.
#'
#' @param power_results Simulation results from a `brms_inla_power*` function with `compute_bayes_factor = TRUE`.
#' @param x_effect Name of effect grid column for x-axis (default: first grid column).
#' @param y_effect Name of effect grid column for y-axis (default: "n").
#' @param facet_by Optional column(s) to facet by.
#' @param n Optional sample size to filter to (NULL means plot all; else show only that n).
#' @param agg_fun Aggregation function (default: mean).
#' @param title,subtitle Optional plot labels.
#' @return ggplot object.
#' @export
plot_bf_heatmap <- function(
    power_results,
    x_effect = NULL,
    y_effect = "n",
    facet_by = NULL,
    n = NULL,
    agg_fun = mean,
    title = NULL,
    subtitle = NULL
) {
  s <- power_results$summary
  # Detect effect grid columns
  effect_names <- setdiff(
    colnames(s),
    c("n", "power_direction", "power_threshold", "power_rope",
      "avg_ci_width", "ci_coverage", "ciw_q05", "ciw_q25", "ciw_q50", "ciw_q75", "ciw_q95",
      "avg_post_prob_direction", "avg_post_prob_threshold", "avg_post_prob_rope",
      "bf_hit_3", "bf_hit_10", "mean_log10_bf", "nsims_ok")
  )
  x_effect <- x_effect %||% effect_names[1]
  y_effect <- y_effect %||% ("n" %in% names(s) ? "n" : effect_names[2] %||% effect_names[1])
  stopifnot(x_effect %in% names(s), y_effect %in% names(s))
  facet_names <- facet_by %||% character(0)
  facet_names <- intersect(facet_names, colnames(s))
  plot_dat <- s
  # If requested, filter for a fixed n
  if (!is.null(n) && "n" %in% names(plot_dat)) plot_dat <- plot_dat[plot_dat$n == n, , drop = FALSE]
  # Aggregate duplicates (if needed)
  group_syms <- c(x_effect, y_effect, facet_names)
  plot_dat <- plot_dat %>%
    dplyr::group_by(!!!rlang::syms(group_syms)) %>%
    dplyr::summarise(val = agg_fun(mean_log10_bf, na.rm = TRUE), .groups = "drop")
  # Build heatmap (tile)
  p <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(x = !!rlang::sym(x_effect), y = !!rlang::sym(y_effect), fill = val)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "C", name = "E[log10 BF10]") +
    ggplot2::labs(
      x = x_effect,
      y = y_effect,
      title = title %||% "Expected Evidence Heatmap",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  # Faceting support
  if (length(facet_names) == 1) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_names)))
  } else if (length(facet_names) > 1) {
    p <- p + ggplot2::facet_grid(as.formula(paste(paste(facet_names, collapse = " + "), "~ .")))
  }
  return(p)
}


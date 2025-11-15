
#' Plot Precision Assurance Curve (Multi-Effect Grid Friendly)
#'
#' Plots the assurance (proportion of runs meeting CI width <= target) vs. a chosen effect grid variable across sample size(s).
#' Supports faceting, effect filtering, and weights.
#'
#' @param power_results List returned by `brms_inla_power*`.
#' @param precision_target Numeric; credible interval width threshold for success.
#' @param x_effect Name of effect grid column for x-axis (default: first grid column).
#' @param facet_by Optional effect grid column(s) for faceting.
#' @param effect_filters Optional named list for filtering rows, e.g. list(treatment=0).
#' @param effect_weights Optional named numeric vector for weights over selected x_effect values.
#' @param title,subtitle Optional plot labels.
#' @return A ggplot object.
#' @export
plot_precision_assurance_curve <- function(
    power_results,
    precision_target,
    x_effect = NULL,
    facet_by = NULL,
    effect_filters = NULL,
    effect_weights = NULL,
    title = NULL,
    subtitle = NULL
) {
  df <- dplyr::filter(power_results$results, ok)
  # Auto-detect effect grid columns
  effect_names <- setdiff(
    colnames(df),
    c("n", "ok", "ci_width", "ci_lower", "ci_upper",
      "post_prob_direction", "post_prob_threshold", "post_prob_rope",
      "bf10", "log10_bf10")
  )
  # Apply user filters (e.g. treatment==0)
  if (!is.null(effect_filters) && length(effect_filters)) {
    for (ef in names(effect_filters)) {
      df <- df[df[[ef]] == effect_filters[[ef]], , drop = FALSE]
    }
  }
  x_effect <- x_effect %||% effect_names[1]
  stopifnot(!is.null(x_effect), x_effect %in% names(df))
  facet_names <- facet_by %||% character(0)
  facet_names <- intersect(facet_names, colnames(df))
  # Group by sample size, x_effect, facets
  tmp <- df %>%
    dplyr::group_by(n, !!rlang::sym(x_effect), !!!rlang::syms(facet_names)) %>%
    dplyr::summarise(
      assurance = mean(ci_width <= precision_target, na.rm = TRUE),
      .groups = "drop"
    )
  # Optional: weighted aggregation over x_effect
  if (!is.null(effect_weights) && is.numeric(effect_weights)) {
    w <- tibble::tibble(val = as.numeric(names(effect_weights)), w = as.numeric(effect_weights))
    names(w)[1] <- x_effect
    tmp <- dplyr::left_join(tmp, w, by = x_effect) %>%
      dplyr::group_by(n, !!!rlang::syms(facet_names)) %>%
      dplyr::summarise(
        assurance = stats::weighted.mean(assurance, w, na.rm = TRUE),
        .groups = "drop"
      )
  }
  p <- ggplot2::ggplot(
    tmp,
    ggplot2::aes(x = !!rlang::sym(x_effect), y = assurance, group = n, color = factor(n))
  ) +
    .geom_line_lw(width = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_viridis_d(name = "Sample size") +
    ggplot2::scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = x_effect,
      y = paste0("Assurance: P(width <=", precision_target, ")"),
      title = title %||% "Precision assurance curve",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  if (length(facet_names) == 1) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_names)))
  } else if (length(facet_names) > 1) {
    p <- p + ggplot2::facet_grid(as.formula(paste(paste(facet_names, collapse = " + "), "~ .")))
  }
  return(p)
}


#' Precision assurance as a function of sample size
#'
#' Plots the proportion of simulations in which the posterior credible
#' interval width is less than or equal to a target, as a function of
#' sample size n. Optionally colours separate curves by an effect-grid
#' variable.
#'
#' This implementation works directly from the per-simulation results
#' (column `ci_width`) and does not rely on the robustness engine.
#'
#' @param power_results Output from a `brms_inla_power*` function, or a
#'   data.frame with at least columns `n` and `ci_width`, plus any
#'   effect-grid columns (e.g. `treatment`, `age_effect`).
#' @param ci_width_target Numeric, target width for the credible interval.
#'   Assurance is defined as Pr(ci_width <= ci_width_target).
#' @param effect_filter Optional named list for filtering effect-grid
#'   columns, e.g. `list(treatment = 0.3)`.
#' @param colour_by Optional name of an effect-grid column to colour
#'   separate curves by. If `NULL`, only n is used.
#' @param title,subtitle Optional plot labels.
#'
#' @return A ggplot object.
#' @export
plot_precision_fan_chart <- function(
    power_results,
    ci_width_target,
    effect_filter = NULL,
    colour_by     = NULL,
    title         = NULL,
    subtitle      = NULL
) {
  # Extract per-simulation data frame
  df <- if (is.list(power_results) && !is.null(power_results$results)) {
    power_results$results
  } else {
    power_results
  }
  stopifnot(is.data.frame(df))
  
  if (!all(c("n", "ci_width") %in% names(df))) {
    stop("Data must contain columns 'n' and 'ci_width'.", call. = FALSE)
  }
  if (!is.numeric(ci_width_target) || length(ci_width_target) != 1L) {
    stop("ci_width_target must be a single numeric value.", call. = FALSE)
  }
  
  # Identify effect columns = everything not obviously sim-level
  non_eff <- c(
    "sim", "n", "ok", "ci_width",
    grep("^(mean_|sd_)", names(df), value = TRUE),
    grep("^post_prob_", names(df), value = TRUE),
    "bf10", "log10_bf10"
  )
  eff_cols <- setdiff(names(df), non_eff)
  
  # Apply optional effect filters
  if (!is.null(effect_filter) && length(effect_filter)) {
    for (nm in names(effect_filter)) {
      if (!nm %in% names(df)) {
        stop("effect_filter refers to unknown column: ", nm, call. = FALSE)
      }
      df <- df[df[[nm]] == effect_filter[[nm]], , drop = FALSE]
    }
  }
  if (nrow(df) == 0L) {
    stop("No rows left after filtering.", call. = FALSE)
  }
  
  # Decide which column (if any) to use for colour
  if (!is.null(colour_by)) {
    if (!colour_by %in% eff_cols) {
      stop("colour_by '", colour_by, "' is not an effect-grid column.", call. = FALSE)
    }
    colour_var <- colour_by
  } else {
    colour_var <- NULL
  }
  
  # Grouping columns: always by n; optionally by colour_var
  group_cols <- "n"
  if (!is.null(colour_var)) {
    group_cols <- c(group_cols, colour_var)
  }
  
  # Aggregate to get precision assurance with Wilson interval
  wilson <- function(h, t, conf = 0.95) {
    if (is.na(h) || is.na(t) || t == 0) {
      return(c(p = NA_real_, lo = NA_real_, hi = NA_real_))
    }
    z <- qnorm(1 - (1 - conf) / 2)
    p <- h / t
    denom  <- 1 + z^2 / t
    center <- (p + z^2 / (2 * t)) / denom
    half   <- (z * sqrt(p * (1 - p) / t + z^2 / (4 * t^2))) / denom
    c(
      p  = p,
      lo = max(0, center - half),
      hi = min(1, center + half)
    )
  }
  
  agg <- df %>%
    dplyr::mutate(hit = as.integer(is.finite(ci_width) & ci_width <= ci_width_target)) %>%
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
  
  # Build the plot
  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(x = n, y = p)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lo, ymax = hi),
      alpha  = 0.15,
      colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point()
  
  if (!is.null(colour_var)) {
    p <- p +
      ggplot2::aes(colour = factor(.data[[colour_var]])) +
      ggplot2::scale_colour_viridis_d(name = colour_var)
  }
  
  p +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x        = "n",
      y        = "Assurance P(width <= target)",
      title    = title %||% "Precision assurance vs sample size",
      subtitle = subtitle %||% paste0("Target CI width: ", ci_width_target)
    ) +
    ggplot2::theme_minimal()
}
#' Add sample-size decision overlay to an assurance contour
#'
#' Overlays the output of `decide_sample_size()` as a step-line
#' on a contour plot (e.g. from `plot_power_contour()`).
#'
#' If the decisions contain multiple effect-grid columns, the
#' overlay is aggregated over all effect columns except `x_effect`
#' by taking the worst-case (maximum) recommended n at each value
#' of `x_effect`.
#'
#' @param p A ggplot object, typically from `plot_power_contour()`.
#' @param decisions A data.frame or tibble returned by
#'   `decide_sample_size()`, containing at least one effect column
#'   and `n_recommended`.
#' @param x_effect Name of the effect column to use on the x-axis.
#'   If `NULL`, the first column that is not `n_recommended`
#'   or `rationale` is used.
#' @param colour Colour for the overlay line (default "red").
#'
#' @return A ggplot object.
#' @export
add_decision_overlay <- function(
    p,
    decisions,
    x_effect = NULL,
    colour   = "red"
) {
  stopifnot(inherits(p, "ggplot"))
  stopifnot(is.data.frame(decisions))
  
  # Determine x_effect if not supplied
  if (is.null(x_effect)) {
    candidates <- setdiff(names(decisions), c("n_recommended", "rationale"))
    if (length(candidates) == 0L) {
      stop("Could not infer x_effect column from decisions.", call. = FALSE)
    }
    x_effect <- candidates[1L]
  }
  if (!x_effect %in% names(decisions)) {
    stop("x_effect column '", x_effect, "' not found in decisions.", call. = FALSE)
  }
  
  # Drop rows with NA recommended n
  df <- decisions[!is.na(decisions$n_recommended), , drop = FALSE]
  if (nrow(df) == 0L) {
    # Nothing meets the targets; return plot unchanged
    return(p)
  }
  
  # Aggregate over any other effect-grid dimensions by worst case (max n)
  other_effects <- setdiff(names(df), c(x_effect, "n_recommended", "rationale"))
  if (length(other_effects) > 0L) {
    df <- df |>
      dplyr::group_by(!!rlang::sym(x_effect)) |>
      dplyr::summarise(
        n_recommended = max(n_recommended, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # Ensure sorted by x for a sensible step path
  df <- df[order(df[[x_effect]]), , drop = FALSE]
  
  x_sym <- rlang::sym(x_effect)
  
  p +
    ggplot2::geom_step(
      data        = df,
      mapping     = ggplot2::aes(x = !!x_sym, y = n_recommended),
      inherit.aes = FALSE,
      colour      = colour,
      linewidth   = 0.8
    )
}
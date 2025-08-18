#' Plot Decision Assurance Curve (Multi-Effect Grid Friendly)
#'
#' Plots the assurance (proportion of simulation runs meeting a
#' posterior probability decision rule) versus an effect grid variable,
#' for a given metric ("direction", "threshold", or "rope") at a fixed
#' decision probability threshold p_star.
#'
#' @param power_results A list returned by `brms_inla_power*`.
#' @param metric Decision metric: "direction", "threshold", or "rope".
#' @param p_star Numeric decision threshold in (0,1).
#' @param x_effect Name of effect grid column for x-axis (default: first grid column).
#' @param facet_by Optional effect grid column(s) to facet by.
#' @param effect_filters Optional named list for filtering rows, e.g. list(treatment=0).
#' @param effect_weights Optional named numeric vector of weights for selected x_effect values.
#' @param title,subtitle Optional plot labels.
#' @return A ggplot object.
#' @export
plot_decision_assurance_curve <- function(
    power_results,
    metric = c("direction","threshold","rope"),
    p_star = 0.95,
    x_effect = NULL,
    facet_by = NULL,
    effect_filters = NULL,
    effect_weights = NULL,
    title = NULL,
    subtitle = NULL
) {
  metric <- match.arg(metric)
  colname <- switch(metric,
                    direction = "post_prob_direction",
                    threshold = "post_prob_threshold",
                    rope = "post_prob_rope"
  )
  df <- dplyr::filter(power_results$results, ok)
  # Detect effect grid columns
  effect_names <- setdiff(
    colnames(df),
    c("n", "ok", "post_prob_direction", "post_prob_threshold", "post_prob_rope",
      "ci_width", "ci_lower", "ci_upper", "bf10", "log10_bf10")
  )
  # Optional: user-restricted grid rows
  if (!is.null(effect_filters) && length(effect_filters)) {
    for (ef in names(effect_filters)) {
      df <- df[df[[ef]] == effect_filters[[ef]], , drop = FALSE]
    }
  }
  # Main x-axis effect (default: first grid column)
  x_effect <- x_effect %||% effect_names[1]
  stopifnot(!is.null(x_effect), x_effect %in% names(df))
  facet_names <- facet_by %||% character(0)
  facet_names <- intersect(facet_names, colnames(df))
  # Group by: sample size + x_effect + facet(s)
  tmp <- df %>%
    dplyr::group_by(n, !!rlang::sym(x_effect), !!!rlang::syms(facet_names)) %>%
    dplyr::summarise(
      assurance = mean(.data[[colname]] >= p_star, na.rm = TRUE),
      .groups = "drop"
    )
  # If requested: weighted means over x_effect
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
  # Plot: lines by sample size and possible faceting
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
      y = paste0("Assurance: P{ post-prob >= ", p_star, " }"),
      title = title %||% paste("Decision assurance (", metric, ")", sep = ""),
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


#' Plot Decision Threshold Contour (Multi-Effect Grid Friendly)
#'
#' Shows assurance as a function of decision threshold p* and one effect grid column, optionally faceted.
#'
#' @param power_results brms_inla_power list (or two-stage, etc.)
#' @param metric Which metric: "direction", "threshold", "rope"
#' @param p_star_grid Numeric vector of decision thresholds (default: 0.5 to 0.99 by 0.01)
#' @param effect_var Name of effect grid column for y-axis (default: first detected grid column)
#' @param facet_by Optional effect grid column(s) to facet by
#' @param effect_value Optional value(s) to filter for effect_var, or named list for multi-filter
#' @param effect_weights Optional weights for aggregation (named by effect_var values)
#' @param title,subtitle Optional plot labels.
#' @return ggplot2 object.
#' @export
plot_decision_threshold_contour <- function(
    power_results,
    metric = c("direction","threshold","rope"),
    p_star_grid = seq(0.5, 0.99, by = 0.01),
    effect_var = NULL,
    facet_by = NULL,
    effect_value = NULL,
    effect_weights = NULL,
    title = NULL,
    subtitle = NULL
) {
  metric <- match.arg(metric)
  colname <- switch(metric,
                    direction = "post_prob_direction",
                    threshold = "post_prob_threshold",
                    rope = "post_prob_rope"
  )
  df <- dplyr::filter(power_results$results, ok)
  effect_names <- setdiff(
    colnames(df),
    c("n", "ok", "post_prob_direction", "post_prob_threshold", "post_prob_rope",
      "ci_width", "ci_lower", "ci_upper", "bf10", "log10_bf10")
  )
  effect_var <- effect_var %||% effect_names[1]
  stopifnot(effect_var %in% names(df))
  # Optionally filter effect_var(s)
  if (!is.null(effect_value)) {
    if (!is.list(effect_value)) effect_value <- setNames(list(effect_value), effect_var)
    for (ef in names(effect_value)) {
      df <- df[df[[ef]] == effect_value[[ef]], , drop = FALSE]
    }
  }
  facet_list <- facet_by %||% character(0)
  facet_names <- intersect(facet_list, colnames(df))
  # Build grid over p* thresholds
  grid <- lapply(p_star_grid, function(ps) {
    tmp <- df %>%
      dplyr::group_by(n, !!rlang::sym(effect_var), !!!rlang::syms(facet_names)) %>%
      dplyr::summarise(
        assurance = mean(.data[[colname]] >= ps, na.rm = TRUE),
        .groups = "drop"
      )
    if (!is.null(effect_weights)) {
      w <- tibble::tibble(val = as.numeric(names(effect_weights)), w = as.numeric(effect_weights))
      names(w)[1] <- effect_var
      tmp <- dplyr::left_join(tmp, w, by = effect_var) %>%
        dplyr::group_by(n, !!!rlang::syms(facet_names)) %>%
        dplyr::summarise(
          assurance = stats::weighted.mean(assurance, w, na.rm = TRUE),
          .groups = "drop"
        )
    }
    dplyr::mutate(tmp, p_star = ps)
  }) %>% dplyr::bind_rows()
  # Contour plot: x=p*, y=effect, z=assurance
  p <- ggplot2::ggplot(grid,
                       ggplot2::aes(x = p_star, y = !!rlang::sym(effect_var), z = assurance)) +
    ggplot2::geom_contour_filled(bins = 12, alpha = 0.9) +
    .add_contour_lines(colour = "white", alpha = 0.3, width = 0.2) +
    ggplot2::scale_fill_viridis_d(name = "Assurance") +
    ggplot2::labs(
      x = "Decision threshold p*",
      y = effect_var,
      title = title %||% paste("Assurance contour for", metric),
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

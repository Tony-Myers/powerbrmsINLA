
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



#' Plot Precision Assurance Fan Chart (Multi-Effect Grid Friendly)
#'
#' Shows assurance (proportion of runs meeting CI width <= target) across sample size(s) and effect grid.
#' Optionally overlays the range (fan/ribbon) across multiple scenarios.
#'
#' @param power_results_list Named list of brms_inla_power results (for fan chart across scenarios) or a single object.
#' @param ci_width_target Numeric; target credible interval width.
#' @param x_effect Name of effect grid column for x-axis (default: first detected grid column).
#' @param facet_by Optional grid column(s) for faceting.
#' @param effect_filters Optional named list for filtering rows, e.g. list(treatment=0).
#' @param effect_weights Optional named numeric vector for averaging over grid values.
#' @param show_individual_scenarios Logical: overlay all scenario curves if TRUE.
#' @param title,subtitle Optional plot labels.
#' @return A ggplot object.
#' @export
plot_precision_fan_chart <- function(
    power_results_list,
    ci_width_target,
    x_effect = NULL,
    facet_by = NULL,
    effect_filters = NULL,
    effect_weights = NULL,
    show_individual_scenarios = FALSE,
    title = NULL,
    subtitle = NULL
) {
  # Accept either a single results object or list of scenarios
  if (!is.list(power_results_list) || !all(sapply(power_results_list, function(x) "results" %in% names(x)))) {
    power_results_list <- list(Default = power_results_list)
  }
  scenario_names <- names(power_results_list)
  if (is.null(scenario_names) || any(scenario_names == "")) {
    scenario_names <- paste0("Scenario ", seq_along(power_results_list))
  }
  per_scen <- lapply(seq_along(power_results_list), function(i) {
    pr <- power_results_list[[i]]
    df <- pr$results %>% dplyr::filter(ok)
    effect_names <- setdiff(
      colnames(df),
      c("n", "ok", "post_prob_direction", "post_prob_threshold", "post_prob_rope", "bf10", "ci_width",
        "ci_lower", "ci_upper", "log10_bf10")
    )
    if (!is.null(effect_filters) && length(effect_filters)) {
      for (ef in names(effect_filters)) {
        df <- df[df[[ef]] == effect_filters[[ef]], , drop = FALSE]
      }
    }
    x_effect_cur <- x_effect %||% effect_names[1]
    facet_names <- facet_by %||% character(0)
    facet_names <- intersect(facet_names, colnames(df))
    tmp <- dplyr::group_by(df, n, !!rlang::sym(x_effect_cur), !!!rlang::syms(facet_names)) %>%
      dplyr::summarise(assurance = mean(ci_width <= ci_width_target, na.rm = TRUE), .groups = "drop")
    if (!is.null(effect_weights)) {
      w <- tibble::tibble(val = as.numeric(names(effect_weights)), w = as.numeric(effect_weights))
      names(w)[1] <- x_effect_cur
      tmp <- dplyr::left_join(tmp, w, by = x_effect_cur) %>%
        dplyr::group_by(n, !!!rlang::syms(facet_names)) %>%
        dplyr::summarise(
          assurance = stats::weighted.mean(assurance, w, na.rm = TRUE),
          .groups = "drop"
        )
    }
    tmp$scenario <- scenario_names[i]
    tmp
  }) %>% dplyr::bind_rows()
  # Ribbon: min-max assurance across scenarios per sample size/grid combo
  group_syms <- unique(c("n", x_effect %||% names(per_scen)[1], intersect(facet_by, colnames(per_scen))))
    dplyr::group_by(!!!rlang::syms(group_syms)) %>%
    dplyr::summarise(
      ymin = min(assurance, na.rm = TRUE),
      ymax = max(assurance, na.rm = TRUE),
      ymid = mean(assurance, na.rm = TRUE),
      .groups = "drop"
    )
  # Plot
  p <- ggplot2::ggplot(ribbon, ggplot2::aes(
    x = !!rlang::sym(x_effect %||% names(per_scen)[1]),
    y = ymid, group = n, color = factor(n))) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), alpha = 0.15) +
    .geom_line_lw(width = 1)
  if (isTRUE(show_individual_scenarios)) {
    p <- p +
      .geom_line_lw(data = per_scen,
                    mapping = ggplot2::aes(y = assurance, color = scenario),
                    width = 0.8)
    if ("scale_colour_viridis_d" %in% getNamespaceExports("ggplot2")) {
      p <- p + ggplot2::scale_colour_viridis_d(name = "Scenario", option = "plasma")
    } else {
      p <- p + ggplot2::scale_colour_manual(values = viridisLite::viridis(length(unique(per_scen$scenario))))
    }
  } else {
    p <- p + ggplot2::scale_color_viridis_d(name = "Sample size")
  }
  p <- p +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = x_effect %||% names(per_scen)[1],
      y = paste0("Assurance P(width <= ", ci_width_target, ")"),
      title = title %||% "Precision assurance (fan chart/ribbon)",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  if (!is.null(facet_by) && length(facet_by) == 1) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_by)))
  } else if (!is.null(facet_by) && length(facet_by) > 1) {
    p <- p + ggplot2::facet_grid(as.formula(paste(paste(facet_by, collapse = " + "), "~ .")))
  }
  return(p)
}

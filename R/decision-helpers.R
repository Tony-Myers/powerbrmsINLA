#' Decide recommended sample size from power/assurance results
#'
#' Returns the smallest per-group sample size that meets user-specified power or
#' assurance targets.  Works with [brms_inla_power()] and related engine
#' outputs.
#'
#' Targets may be supplied as direct arguments or via a `targets` list; direct
#' arguments take precedence.
#'
#' ## Modes of operation
#'
#' **Assurance mode** (recommended): supply `prior_weights` with a design prior
#' over the effect grid.  The function calls [compute_assurance()] internally
#' and returns the smallest sample size where *unconditional* Bayesian
#' assurance reaches `target` for each requested metric.
#'
#' **Conditional mode** (backward compatible): when `prior_weights = NULL`, the
#' function selects the smallest n at which the per-cell conditional power meets
#' the specified target for each effect-size value.  If the summary contains
#' sampled variance columns (e.g., from distributional `error_sd`), they are
#' averaged out before selecting n.
#'
#' @param x A list with `$summary` (engine output) or a plain `data.frame`
#'   summary.
#' @param direction Numeric in \eqn{[0,1]}: required power for
#'   `power_direction` (optional).
#' @param threshold Numeric in \eqn{[0,1]}: required power for
#'   `power_threshold` (optional).
#' @param rope_in Numeric in \eqn{[0,1]}: maximum allowed Pr(in ROPE)
#'   (optional).
#' @param bf10 Numeric Bayes-factor cutoff (e.g., 10).  If provided the
#'   function looks for a column `bf_hit_<bf10>`; if absent it falls back to
#'   per-simulation bf10 in `x$results`.
#' @param bf_prop_min Numeric in \eqn{[0,1]}: minimum proportion of simulations
#'   that must achieve BF ≥ `bf10` (default 0; conditional mode only).
#' @param targets Optional named list alternative to the direct arguments.
#'   Ignored when any direct argument is non-`NULL`.
#' @param prior_weights A design prior for assurance mode.  Accepts the same
#'   formats as [compute_assurance()]: a named numeric vector of weights over
#'   effect-grid values, or a distribution list such as
#'   `list(dist = "normal", mean = 0.5, sd = 0.2)`.  When `NULL` (default),
#'   the function operates in conditional mode.
#' @param target Numeric in \eqn{[0,1]}: assurance level to achieve (default
#'   `0.80`; assurance mode only).
#'
#' @return An object of class `"powerbrmsINLA_sample_size"` (which inherits from
#'   `data.frame`) with a [print.powerbrmsINLA_sample_size()] method.
#'
#'   **Assurance mode** columns:
#'   \describe{
#'     \item{`metric`}{Requested metric name (`"direction"`, `"threshold"`,
#'       `"rope_in"`, `"bf10"`).}
#'     \item{`target`}{The assurance target supplied.}
#'     \item{`n_recommended`}{Smallest per-group sample size achieving the
#'       target, or `NA` if none in the grid qualifies.}
#'     \item{`assurance_achieved`}{Assurance value at the recommended n.}
#'     \item{`prior_description`}{Plain-text description of the design prior.}
#'   }
#'
#'   **Conditional mode** columns:
#'   \describe{
#'     \item{`<effect column(s)>`}{One column per effect variable.}
#'     \item{`n_recommended`}{Smallest per-group sample size meeting all
#'       targets, or `NA`.}
#'     \item{`cond_power_*`}{Conditional power at the recommended n for each
#'       requested metric.}
#'   }
#'
#' @seealso [compute_assurance()], [assurance_prior_weights()],
#'   [beta_weights_on_grid()]
#'
#' @export
#'
#' @examples
#' # Build a small synthetic power_result without running INLA
#' syn_summary <- data.frame(
#'   n               = rep(c(50, 100, 200), each = 3),
#'   treatment       = rep(c(0.2, 0.5, 0.8), 3),
#'   power_direction = c(0.40, 0.65, 0.85,
#'                       0.60, 0.82, 0.95,
#'                       0.72, 0.90, 0.98),
#'   power_threshold = c(0.30, 0.55, 0.75,
#'                       0.50, 0.72, 0.88,
#'                       0.62, 0.80, 0.92),
#'   stringsAsFactors = FALSE
#' )
#' syn_result <- list(
#'   summary  = syn_summary,
#'   settings = list(effect_name = "treatment")
#' )
#'
#' # --- Assurance mode ---
#' w <- assurance_prior_weights(c(0.2, 0.5, 0.8), dist = "normal",
#'                               mean = 0.5, sd = 0.2)
#' rec_assurance <- decide_sample_size(
#'   syn_result,
#'   direction     = 0.80,
#'   threshold     = 0.75,
#'   prior_weights = w,
#'   target        = 0.80
#' )
#' print(rec_assurance)
#'
#' # --- Conditional mode (backward compatible) ---
#' rec_cond <- decide_sample_size(syn_result, direction = 0.80)
#' print(rec_cond)
decide_sample_size <- function(
    x,
    direction     = NULL,
    threshold     = NULL,
    rope_in       = NULL,
    bf10          = NULL,
    bf_prop_min   = 0,
    targets       = NULL,
    prior_weights = NULL,
    target        = 0.80
) {
  # Merge explicit targets with list; direct arguments take precedence
  if (is.null(targets)) targets <- list()
  if (!is.null(direction)) targets$direction <- direction
  if (!is.null(threshold)) targets$threshold <- threshold
  if (!is.null(rope_in))   targets$rope_in   <- rope_in
  if (!is.null(bf10))      targets$bf10      <- bf10

  if (!is.null(prior_weights)) {
    return(.decide_assurance_mode(x, targets, prior_weights, target))
  }

  .decide_conditional_mode(x, targets, bf_prop_min)
}


# ---------------------------------------------------------------------------
# Assurance mode internals
# ---------------------------------------------------------------------------

.decide_assurance_mode <- function(x, targets, prior_weights, target) {
  # Map from decide_sample_size target names to compute_assurance metric names
  metric_map <- c(
    direction = "direction",
    threshold = "threshold",
    rope_in   = "rope",
    bf10      = "bf"
  )

  requested <- intersect(names(targets), names(metric_map))

  if (length(requested) == 0L) {
    stop(
      "No metric targets supplied. Provide at least one of: ",
      "direction, threshold, rope_in, bf10.",
      call. = FALSE
    )
  }

  prior_desc <- .format_prior_description(prior_weights)

  rows <- lapply(requested, function(tgt_name) {
    ca_metric <- metric_map[[tgt_name]]

    assur_obj <- tryCatch(
      compute_assurance(x, prior_weights = prior_weights, metric = ca_metric),
      error = function(e) {
        message(
          "decide_sample_size(): cannot compute assurance for '",
          tgt_name, "': ", conditionMessage(e)
        )
        NULL
      }
    )

    if (is.null(assur_obj)) {
      return(data.frame(
        metric             = tgt_name,
        target             = target,
        n_recommended      = NA_integer_,
        assurance_achieved = NA_real_,
        prior_description  = prior_desc,
        stringsAsFactors   = FALSE
      ))
    }

    df <- assur_obj$assurance
    df <- df[order(df$sample_size), , drop = FALSE]
    ok <- is.finite(df$assurance) & (df$assurance >= target)

    if (any(ok)) {
      n_rec <- min(df$sample_size[ok])
      a_ach <- df$assurance[df$sample_size == n_rec][1L]
    } else {
      message(
        "decide_sample_size(): no sample size in the grid achieves ",
        round(target * 100), "% assurance for metric '", tgt_name,
        "'. Returning NA."
      )
      n_rec <- NA_integer_
      a_ach <- NA_real_
    }

    data.frame(
      metric             = tgt_name,
      target             = target,
      n_recommended      = as.integer(n_rec),
      assurance_achieved = a_ach,
      prior_description  = prior_desc,
      stringsAsFactors   = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("powerbrmsINLA_sample_size", "data.frame")
  attr(out, "mode") <- "assurance"
  out
}


# ---------------------------------------------------------------------------
# Conditional mode internals
# ---------------------------------------------------------------------------

.decide_conditional_mode <- function(x, targets, bf_prop_min) {
  s <- .get_summary_df(x)
  if (!"n" %in% names(s)) {
    stop("Summary data must contain a column 'n'.", call. = FALSE)
  }
  if (nrow(s) == 0L) {
    stop("Summary data has 0 rows; cannot decide sample size.", call. = FALSE)
  }
  s <- as.data.frame(s)
  s <- s[order(s$n), , drop = FALSE]

  # Identify effect columns; exclude sampled-SD columns so they are not
  # treated as design factors
  non_eff <- unique(c(
    "n",
    "power_direction", "power_threshold", "power_rope",
    "avg_ci_width", "ci_coverage", "ci_width", "ci_lower", "ci_upper",
    "bf10_mean", "bf10_sd", "bf_hit_3", "bf_hit_10",
    grep("^bf_hit_",  names(s), value = TRUE),
    grep("^mean_log", names(s), value = TRUE),
    "bf_median", "bf_min", "bf_max", "mean_log10_bf", "nsims_ok",
    grep("^sampled_", names(s), value = TRUE)
  ))
  eff_cols <- setdiff(names(s), non_eff)

  # Aggregate across sampled-SD draws if those columns are present
  sampled_cols <- grep("^sampled_", names(s), value = TRUE)
  if (length(sampled_cols) > 0L && length(eff_cols) > 0L) {
    pwr_agg <- intersect(
      c("power_direction", "power_threshold", "power_rope",
        grep("^bf_hit_", names(s), value = TRUE)),
      names(s)
    )
    if (length(pwr_agg) > 0L) {
      group_keys <- c("n", eff_cols)
      s <- as.data.frame(
        dplyr::summarise(
          dplyr::group_by(s, dplyr::across(dplyr::all_of(group_keys))),
          dplyr::across(dplyr::all_of(pwr_agg), \(v) mean(v, na.rm = TRUE)),
          .groups = "drop"
        )
      )
      s <- s[order(s$n), , drop = FALSE]
    }
  }

  ## --- BF proportion ---------------------------------------------------------
  s$bf_prop <- NA_real_

  if (!is.null(targets$bf10)) {
    exact_col <- paste0("bf_hit_", as.integer(targets$bf10))

    if (exact_col %in% names(s)) {
      s$bf_prop <- s[[exact_col]]
    } else {
      r <- NULL
      if (is.list(x) && !is.null(x$results)) r <- x$results

      if (!is.null(r) && all(c("n", "bf10") %in% names(r))) {
        r <- as.data.frame(r)
        eff_cols_r   <- intersect(eff_cols, names(r))
        group_cols_r <- c("n", eff_cols_r)

        bf_df <- r |>
          dplyr::group_by(dplyr::across(dplyr::all_of(group_cols_r))) |>
          dplyr::summarise(
            bf_prop = mean(bf10 >= targets$bf10, na.rm = TRUE),
            .groups = "drop"
          )

        s <- dplyr::left_join(
          s, bf_df,
          by    = c("n", eff_cols_r),
          suffix = c("", "_from_results")
        )

        if ("bf_prop_from_results" %in% names(s)) {
          s$bf_prop               <- s$bf_prop_from_results
          s$bf_prop_from_results  <- NULL
          message(
            "decide_sample_size(): recomputed Pr(BF10 >= ",
            targets$bf10, ") from per-simulation bf10 values."
          )
        } else {
          message(
            "decide_sample_size(): could not align per-simulation BF values ",
            "with summary; BF target will be ignored."
          )
        }
      } else {
        message(
          "decide_sample_size(): BF target requested but no matching bf_hit_* ",
          "column and no usable per-simulation bf10 in x$results; ",
          "BF target will be ignored."
        )
      }
    }
  }

  ## --- Vectorised criteria ---------------------------------------------------
  ok_dir <- rep(TRUE, nrow(s))
  if (!is.null(targets$direction) && "power_direction" %in% names(s)) {
    ok_dir <- is.finite(s$power_direction) & (s$power_direction >= targets$direction)
  }

  ok_thr <- rep(TRUE, nrow(s))
  if (!is.null(targets$threshold) && "power_threshold" %in% names(s)) {
    ok_thr <- is.finite(s$power_threshold) & (s$power_threshold >= targets$threshold)
  }

  # power_rope = Pr(outside ROPE); rope_in target is Pr(in ROPE) <= rope_in
  ok_rope <- rep(TRUE, nrow(s))
  if (!is.null(targets$rope_in) && "power_rope" %in% names(s)) {
    p_out  <- s$power_rope
    p_in   <- ifelse(is.finite(p_out), 1 - p_out, NA_real_)
    ok_rope <- is.finite(p_in) & (p_in <= targets$rope_in)
  }

  ok_bf <- rep(TRUE, nrow(s))
  if (!is.null(targets$bf10) && any(is.finite(s$bf_prop))) {
    ok_bf <- is.finite(s$bf_prop) & (s$bf_prop >= bf_prop_min)
  }

  ok_all <- ok_dir & ok_thr & ok_rope & ok_bf

  ## --- Group by effect grid and pick minimal n ------------------------------
  if (length(eff_cols) > 0L) {
    s$._group <- interaction(s[eff_cols], drop = TRUE, lex.order = TRUE)
  } else {
    s$._group <- factor("all")
    eff_cols  <- character(0L)
  }
  s$._ok_all <- ok_all

  groups <- levels(droplevels(s$._group))

  rows <- lapply(groups, function(grp) {
    sg <- s[s$._group == grp, , drop = FALSE]

    if (any(sg$._ok_all, na.rm = TRUE)) {
      n_rec    <- min(sg$n[sg$._ok_all], na.rm = TRUE)
      sg_rec   <- sg[sg$n == n_rec, , drop = FALSE][1L, , drop = FALSE]
      no_match <- FALSE
    } else {
      n_rec    <- NA_integer_
      sg_rec   <- sg[1L, , drop = FALSE]
      no_match <- TRUE
    }

    row <- data.frame(n_recommended = as.integer(n_rec), stringsAsFactors = FALSE)

    for (ec in eff_cols) {
      row[[ec]] <- sg[[ec]][1L]
    }

    # Append conditional power at the recommended n (or NA when not met)
    .add_cond_pwr <- function(col_out, col_src, transform = identity) {
      val <- if (!no_match && col_src %in% names(sg_rec)) {
        transform(sg_rec[[col_src]])
      } else {
        NA_real_
      }
      row[[col_out]] <<- val
    }

    if (!is.null(targets$direction) && "power_direction" %in% names(s)) {
      .add_cond_pwr("cond_power_direction", "power_direction")
    }
    if (!is.null(targets$threshold) && "power_threshold" %in% names(s)) {
      .add_cond_pwr("cond_power_threshold", "power_threshold")
    }
    if (!is.null(targets$rope_in) && "power_rope" %in% names(s)) {
      # Report Pr(in ROPE) = 1 - power_rope for interpretability
      .add_cond_pwr("cond_power_rope", "power_rope", transform = function(v) 1 - v)
    }

    row
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  # Column order: effect columns, then n_recommended, then power columns
  pwr_cols  <- intersect(
    c("cond_power_direction", "cond_power_threshold", "cond_power_rope"),
    names(out)
  )
  col_order <- c(eff_cols, "n_recommended", pwr_cols)
  out       <- out[, intersect(col_order, names(out)), drop = FALSE]

  n_na <- sum(is.na(out$n_recommended))
  if (n_na > 0L) {
    message(
      "decide_sample_size(): no sample size met all targets for ",
      n_na, " effect-size group(s). Returning NA."
    )
  }

  class(out) <- c("powerbrmsINLA_sample_size", "data.frame")
  attr(out, "mode") <- "conditional"
  out
}


# ---------------------------------------------------------------------------
# Shared helper: format prior description
# ---------------------------------------------------------------------------

.format_prior_description <- function(prior_weights) {
  if (is.list(prior_weights) && !is.null(prior_weights$dist)) {
    ps        <- prior_weights
    param_str <- paste(
      vapply(setdiff(names(ps), "dist"), function(nm)
        paste0(nm, " = ", ps[[nm]]), character(1L)),
      collapse = ", "
    )
    paste0(ps$dist, "(", param_str, ")")
  } else if (is.numeric(prior_weights)) {
    n_vals <- length(prior_weights)
    nms    <- names(prior_weights)
    if (!is.null(nms) && n_vals <= 6L) {
      paste0(
        "user weights [",
        paste(paste0(nms, "=", round(prior_weights, 3L)), collapse = ", "),
        "]"
      )
    } else {
      paste0("user-supplied weights over ", n_vals, " effect value(s)")
    }
  } else {
    "unknown prior"
  }
}


# ---------------------------------------------------------------------------
# Print method
# ---------------------------------------------------------------------------

#' Print method for powerbrmsINLA_sample_size objects
#'
#' Formats recommendations from [decide_sample_size()] in a human-readable
#' way.  In assurance mode each row is rendered as a plain-English sentence.
#' In conditional mode a formatted table is printed.
#'
#' @param x An object of class `"powerbrmsINLA_sample_size"` returned by
#'   [decide_sample_size()].
#' @param digits Number of decimal places for power/assurance values
#'   (default 4).
#' @param ... Unused; present for S3 compatibility.
#'
#' @return `x`, invisibly.
#' @export
print.powerbrmsINLA_sample_size <- function(x, digits = 4L, ...) {
  mode <- attr(x, "mode") %||% "conditional"

  if (identical(mode, "assurance")) {
    for (i in seq_len(nrow(x))) {
      row <- x[i, ]
      pct <- round(row$target * 100)
      if (is.na(row$n_recommended)) {
        cat(sprintf(
          "No sample size in the grid achieves %d%% assurance using the %s criterion.\n",
          pct, row$metric
        ))
      } else {
        cat(sprintf(
          paste0(
            "Based on %s, the recommended sample size for %d%% assurance ",
            "using the %s criterion is %d per group (assurance = %s).\n"
          ),
          row$prior_description,
          pct,
          row$metric,
          row$n_recommended,
          format(round(row$assurance_achieved, digits), nsmall = digits)
        ))
      }
    }
  } else {
    cat("Conditional power-based sample size recommendations\n")
    cat("===================================================\n")
    df <- x
    class(df) <- "data.frame"
    pwr_cols <- grep("^cond_power_", names(df), value = TRUE)
    for (nm in pwr_cols) {
      df[[nm]] <- round(df[[nm]], digits)
    }
    print(df, row.names = FALSE)
  }

  invisible(x)
}


# ---------------------------------------------------------------------------
# Plot helper (unchanged)
# ---------------------------------------------------------------------------

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

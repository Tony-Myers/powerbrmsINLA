# sequential-analysis.R
# User-facing sequential Bayesian analysis of REAL accumulating data:
#   sequential_design()      - prespecify the design (with fingerprint)
#   sequential_analysis()    - analyse accumulated data at an interim look
#   plot_sequential_monitor()- trajectory plots with stopping boundaries
#
# Companion to brms_inla_sequential_trial(), which simulates the operating
# characteristics of the same design object before data collection.
#
# New in 1.3.0.


# ---------------------------------------------------------------------------
# Design constructor
# ---------------------------------------------------------------------------

#' Prespecify a sequential Bayesian analysis design
#'
#' @description
#' Creates a frozen specification of a sequential Bayesian analysis
#' \emph{before} data collection: the model, the analysis priors, the
#' monitored effect, the decision metric, the stopping thresholds, and the
#' planned look schedule. The returned object is used both to simulate the
#' design's operating characteristics ([brms_inla_sequential_trial()]) and to
#' monitor the real study as data accumulate ([sequential_analysis()]).
#'
#' The object carries an MD5 fingerprint of all decision-relevant fields.
#' Quoting this fingerprint in a preregistration or protocol allows readers
#' and reviewers to verify that the stopping rules were not altered after
#' data collection began: any change to the model, priors, metric,
#' thresholds, or look schedule produces a different fingerprint.
#'
#' @details
#' Because Bayesian posterior quantities obey the likelihood principle, the
#' interpretation of the posterior at the final look does not depend on the
#' stopping rule. The \emph{frequency} properties of the design (false-go
#' rate, expected sample size, exaggeration of estimates at early stops) do.
#' The recommended workflow is therefore: (1) create the design; (2) estimate
#' its operating characteristics by simulation with
#' [brms_inla_sequential_trial()]; (3) preregister the design, including its
#' fingerprint; (4) monitor the study with [sequential_analysis()]; and
#' (5) report the final posterior alongside the simulated operating
#' characteristics.
#'
#' @param formula Model formula (brms syntax; random effects supported).
#' @param family Response family, e.g. `gaussian()`, `binomial()`.
#' @param priors Optional brms prior specification (translated to INLA via
#'   the package's prior bridge; see [brms_inla_power()]).
#' @param effect_name Single character: the monitored fixed effect.
#' @param metric `"direction"`, `"threshold"`, or `"rope"`.
#' @param alternative `"greater"` or `"less"`: the prespecified direction for
#'   the direction/threshold metrics. Ignored for `"rope"`.
#' @param looks Strictly increasing vector of cumulative sample sizes (total
#'   observations) at which interim analyses are planned. The final element
#'   is the maximum sample size.
#' @param prob_success Scalar or `length(looks)` vector of stop-for-success
#'   posterior probability thresholds in (0.5, 1). Per-look vectors allow
#'   stricter early thresholds (O'Brien-Fleming-like behaviour).
#' @param prob_futility Scalar or `length(looks)` vector of stop-for-futility
#'   thresholds in (0, 0.5), or `NULL` (default) for no futility stopping.
#'   Ignored with a warning for `metric = "rope"`.
#' @param effect_threshold Threshold for `metric = "threshold"` (model
#'   coefficient scale).
#' @param rope_bounds Length-2 increasing numeric vector for
#'   `metric = "rope"`.
#' @param credible_level Credible level for reported intervals (default 0.95).
#' @param label Optional character label for the study (used in printing).
#'
#' @return A list of class `"powerbrmsINLA_seq_design"` whose `fingerprint`
#'   element contains the MD5 hash of the decision-relevant fields.
#'
#' @examples
#' design <- sequential_design(
#'   formula          = strength ~ group,
#'   effect_name      = "group",
#'   metric           = "threshold",
#'   effect_threshold = 0.2,
#'   looks            = c(40, 80, 120, 160, 200),
#'   prob_success     = 0.95,
#'   prob_futility    = 0.05,
#'   label            = "Periodisation RCT"
#' )
#' print(design)
#'
#' @seealso [sequential_analysis()], [brms_inla_sequential_trial()],
#'   [plot_sequential_monitor()]
#' @export
sequential_design <- function(
    formula,
    family = gaussian(),
    priors = NULL,
    effect_name,
    metric = c("direction", "threshold", "rope"),
    alternative = c("greater", "less"),
    looks,
    prob_success = 0.95,
    prob_futility = NULL,
    effect_threshold = 0,
    rope_bounds = NULL,
    credible_level = 0.95,
    label = NULL
) {
  metric      <- match.arg(metric)
  alternative <- match.arg(alternative)

  stopifnot(inherits(formula, "formula"))
  stopifnot(is.character(effect_name))
  if (length(effect_name) != 1L) {
    stop("`effect_name` must have length 1: sequential analysis monitors a ",
         "single primary effect.", call. = FALSE)
  }

  if (!is.numeric(looks) || length(looks) < 1L || any(!is.finite(looks)) ||
      any(looks <= 0) || is.unsorted(looks, strictly = TRUE)) {
    stop("`looks` must be a strictly increasing vector of positive sample ",
         "sizes.", call. = FALSE)
  }
  looks   <- as.integer(round(looks))
  n_looks <- length(looks)

  if (!is.numeric(prob_success) ||
      !(length(prob_success) %in% c(1L, n_looks)) ||
      any(prob_success <= 0.5) || any(prob_success >= 1)) {
    stop("`prob_success` must be a scalar or length(looks) vector with ",
         "values in (0.5, 1).", call. = FALSE)
  }
  prob_success <- rep_len(prob_success, n_looks)

  if (!is.null(prob_futility)) {
    if (identical(metric, "rope")) {
      warning("`prob_futility` is ignored for metric = 'rope'; equivalence ",
              "stopping uses `prob_success` applied to P(inside ROPE).",
              call. = FALSE)
      prob_futility <- NULL
    } else {
      if (!is.numeric(prob_futility) ||
          !(length(prob_futility) %in% c(1L, n_looks)) ||
          any(prob_futility <= 0) || any(prob_futility >= 0.5)) {
        stop("`prob_futility` must be NULL, or a scalar / length(looks) ",
             "vector with values in (0, 0.5).", call. = FALSE)
      }
      prob_futility <- rep_len(prob_futility, n_looks)
    }
  }

  if (identical(metric, "rope")) {
    if (is.null(rope_bounds) || length(rope_bounds) != 2L ||
        !is.numeric(rope_bounds) || rope_bounds[1L] >= rope_bounds[2L]) {
      stop("metric = 'rope' requires `rope_bounds` as a length-2 increasing ",
           "numeric vector, e.g. c(-0.2, 0.2).", call. = FALSE)
    }
  }

  stopifnot(is.numeric(credible_level), length(credible_level) == 1L,
            credible_level > 0, credible_level < 1)

  fam_map <- .to_inla_family(family)

  design <- list(
    formula          = formula,
    family           = family,
    inla_family      = fam_map$inla,
    priors           = priors,
    effect_name      = effect_name,
    metric           = metric,
    alternative      = alternative,
    looks            = looks,
    n_looks          = n_looks,
    n_max            = looks[n_looks],
    prob_success     = prob_success,
    prob_futility    = prob_futility,
    effect_threshold = effect_threshold,
    rope_bounds      = rope_bounds,
    credible_level   = credible_level,
    label            = label %||% "Sequential Bayesian analysis",
    created          = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  design$fingerprint <- .seq_design_fingerprint(design)
  class(design) <- c("powerbrmsINLA_seq_design", "list")
  design
}


#' MD5 fingerprint of a sequential design's decision-relevant fields
#'
#' @param design A `powerbrmsINLA_seq_design` list (class not yet required).
#' @return Character scalar MD5 hash.
#' @keywords internal
.seq_design_fingerprint <- function(design) {
  prior_txt <- if (is.null(design$priors)) {
    "default"
  } else {
    pr_df <- tryCatch(as.data.frame(design$priors), error = function(e) NULL)
    if (is.null(pr_df)) "unparseable" else
      paste(apply(pr_df, 1L, paste, collapse = "|"), collapse = ";")
  }
  canon <- paste(
    paste(deparse(design$formula), collapse = " "),
    design$inla_family,
    prior_txt,
    design$effect_name,
    design$metric,
    design$alternative,
    paste(design$looks, collapse = ","),
    paste(format(design$prob_success, digits = 15), collapse = ","),
    if (is.null(design$prob_futility)) "no-futility" else
      paste(format(design$prob_futility, digits = 15), collapse = ","),
    format(design$effect_threshold, digits = 15),
    if (is.null(design$rope_bounds)) "no-rope" else
      paste(format(design$rope_bounds, digits = 15), collapse = ","),
    format(design$credible_level, digits = 15),
    sep = "\n"
  )
  tf <- tempfile(fileext = ".txt")
  on.exit(unlink(tf), add = TRUE)
  writeLines(canon, tf)
  unname(tools::md5sum(tf))
}


#' Print method for sequential design objects
#'
#' @param x A `"powerbrmsINLA_seq_design"` object.
#' @param ... Unused; present for S3 compatibility.
#' @return `x`, invisibly.
#' @export
print.powerbrmsINLA_seq_design <- function(x, ...) {
  cat("Prespecified sequential Bayesian design (powerbrmsINLA)\n")
  cat("========================================================\n")
  cat("Study       :", x$label, "\n")
  cat("Model       :", paste(deparse(x$formula), collapse = " "),
      "[", x$inla_family, "]\n")
  cat("Monitored   :", x$effect_name, "\n")
  cat("Metric      :", x$metric,
      if (!identical(x$metric, "rope"))
        paste0("(alternative: ", x$alternative, ")") else
          paste0("(ROPE: [", x$rope_bounds[1], ", ", x$rope_bounds[2], "])"),
      "\n")
  if (identical(x$metric, "threshold")) {
    cat("Threshold   : effect", if (identical(x$alternative, "greater")) ">" else "<",
        x$effect_threshold, "\n")
  }
  cat("Looks (n)   :", paste(x$looks, collapse = ", "), "\n")
  cat("Stop rules  : success >=",
      paste(format(x$prob_success), collapse = "/"),
      if (!is.null(x$prob_futility))
        paste("; futility <=", paste(format(x$prob_futility), collapse = "/"))
      else "; no futility stopping", "\n")
  cat("Priors      :", if (is.null(x$priors)) "INLA defaults"
      else "brms specification (translated via prior bridge)", "\n")
  cat("Created     :", x$created, "\n")
  cat("Fingerprint :", x$fingerprint, "\n\n")
  cat("Quote the fingerprint in your preregistration/protocol. Estimate the\n")
  cat("design's operating characteristics with brms_inla_sequential_trial()\n")
  cat("before collecting data; monitor the study with sequential_analysis().\n")
  invisible(x)
}


# ---------------------------------------------------------------------------
# Sequential analysis of real data
# ---------------------------------------------------------------------------

#' Analyse accumulated real data at a sequential interim look
#'
#' @description
#' Fits the prespecified model to the data accumulated so far, computes the
#' monitored posterior probability, applies the prespecified stopping
#' threshold for the current look, and records the result in an auditable
#' monitor object. Call it once per interim look, passing the design on the
#' first call and the returned monitor object on subsequent calls.
#'
#' The monitor refuses to continue after a stop decision has been recorded
#' unless `allow_after_stop = TRUE`, in which case the continuation is
#' logged as a protocol deviation. Analyses at sample sizes that do not
#' match the planned schedule are permitted but flagged and recorded.
#'
#' @param x A `"powerbrmsINLA_seq_design"` (first look) or
#'   `"powerbrmsINLA_seq_monitor"` (subsequent looks).
#' @param data Data frame containing all data accumulated so far (not just
#'   the new batch). Must contain the variables in the design formula.
#' @param Ntrials_col Optional name of a column with binomial trial counts
#'   (binomial/beta-binomial families). Defaults to a column named
#'   `.Ntrials` if present, else 1 (Bernoulli).
#' @param allow_after_stop Logical; permit analysis after a recorded stop
#'   decision (logged as a deviation). Default `FALSE`.
#' @param note Optional character note recorded with this look (e.g. reason
#'   for an off-schedule analysis).
#' @param inla_num_threads INLA threads ("outer:inner"); auto-detected if
#'   `NULL`.
#'
#' @return A list of class `"powerbrmsINLA_seq_monitor"` containing the
#'   design, its fingerprint, a `looks` data frame (one row per analysis:
#'   timestamp, n, planned n, estimate, credible interval, monitored
#'   probability, thresholds, decision, deviations), and the current
#'   `status` (`"ongoing"` or a stop decision).
#'
#' @examples
#' \dontrun{
#' design <- sequential_design(
#'   formula = y ~ treatment, effect_name = "treatment",
#'   metric = "direction", looks = c(40, 80, 120),
#'   prob_success = 0.95, prob_futility = 0.05
#' )
#'
#' # After the first 40 observations have been collected:
#' mon <- sequential_analysis(design, data_so_far)
#' print(mon)
#'
#' # After 80:
#' mon <- sequential_analysis(mon, data_so_far)
#' plot_sequential_monitor(mon)
#' }
#'
#' @seealso [sequential_design()], [plot_sequential_monitor()],
#'   [brms_inla_sequential_trial()]
#' @export
sequential_analysis <- function(
    x,
    data,
    Ntrials_col = NULL,
    allow_after_stop = FALSE,
    note = NULL,
    inla_num_threads = NULL
) {
  # ---- Resolve design / monitor -------------------------------------------
  if (inherits(x, "powerbrmsINLA_seq_design")) {
    design  <- x
    monitor <- NULL
  } else if (inherits(x, "powerbrmsINLA_seq_monitor")) {
    monitor <- x
    design  <- x$design
  } else {
    stop("`x` must be a sequential_design() object (first look) or the ",
         "monitor returned by a previous sequential_analysis() call.",
         call. = FALSE)
  }

  # Integrity: the design must not have changed since the monitor was created
  current_fp <- .seq_design_fingerprint(design)
  if (!identical(current_fp, design$fingerprint)) {
    stop("Design fingerprint mismatch: the design object has been modified ",
         "since it was created. Recreate it with sequential_design(); do ",
         "not edit decision rules mid-study.", call. = FALSE)
  }

  stopifnot(is.data.frame(data))
  needed <- setdiff(all.vars(design$formula), c(".Ntrials", ".E"))
  missing_vars <- setdiff(needed, names(data))
  if (length(missing_vars) > 0L) {
    stop("`data` is missing variables required by the design formula: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  # ---- Stop-state handling --------------------------------------------------
  deviations <- character(0L)
  if (!is.null(monitor) && !identical(monitor$status, "ongoing")) {
    if (!isTRUE(allow_after_stop)) {
      stop("This monitor recorded the decision '", monitor$status,
           "' at n = ", monitor$looks$n[nrow(monitor$looks)],
           ". Analysing further data invalidates the prespecified design. ",
           "If you must continue, set allow_after_stop = TRUE; the ",
           "continuation will be recorded as a protocol deviation.",
           call. = FALSE)
    }
    deviations <- c(deviations, "analysis_after_stop")
  }

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required for sequential_analysis(). ",
         "See https://www.r-inla.org for installation instructions.",
         call. = FALSE)
  }

  # ---- Look bookkeeping ------------------------------------------------------
  n <- nrow(data)
  k <- if (is.null(monitor)) 1L else nrow(monitor$looks) + 1L
  k_sched   <- min(k, design$n_looks)
  planned_n <- if (k <= design$n_looks) design$looks[k] else NA_integer_

  if (!is.null(monitor)) {
    last_n <- monitor$looks$n[nrow(monitor$looks)]
    if (n <= last_n) {
      stop("`data` has ", n, " rows but the previous look analysed ",
           last_n, ". Pass ALL accumulated data (previous plus new).",
           call. = FALSE)
    }
  }
  if (!is.na(planned_n) && n != planned_n) {
    deviations <- c(deviations, sprintf("off_schedule(n=%d,planned=%d)", n, planned_n))
  }
  if (k > design$n_looks) {
    deviations <- c(deviations, "unplanned_extra_look")
  }

  ps <- design$prob_success[k_sched]
  pf <- if (is.null(design$prob_futility)) NA_real_ else design$prob_futility[k_sched]

  # ---- Fit with INLA --------------------------------------------------------
  if (is.null(inla_num_threads)) {
    n_cores <- parallel::detectCores()
    if (!is.numeric(n_cores) || length(n_cores) != 1L || !is.finite(n_cores)) n_cores <- 1L
    inla_num_threads <- if (n_cores >= 4) "4:1" else if (n_cores >= 2) "2:1" else "1:1"
  }

  prior_map <- .map_brms_priors_to_inla(
    design$priors,
    family_control_supplied = FALSE,
    inla_family = design$inla_family
  )
  tf_alt <- .brms_to_inla_formula2(design$formula, hyper_by_re = prior_map$hyper_by_re)

  dat <- data
  if (length(tf_alt$re_specs) > 0L) {
    for (re in tf_alt$re_specs) {
      gid <- as.integer(as.factor(dat[[re$group]]))
      if (isTRUE(re$has_intercept) && is.null(dat[[re$id_intercept]]))
        dat[[re$id_intercept]] <- gid
      if (!is.null(re$slope) && is.null(dat[[re$id_slope]]))
        dat[[re$id_slope]] <- gid
    }
  }

  inla_args <- list(
    formula = tf_alt$inla_formula,
    data = dat,
    family = design$inla_family,
    control.fixed = prior_map$control_fixed %||% list(),
    control.predictor = list(link = 1),
    control.family = prior_map$control_family %||% list(),
    num.threads = inla_num_threads,
    verbose = FALSE
  )
  if (design$inla_family %in% c("binomial", "betabinomial")) {
    nt_col <- Ntrials_col %||% if (".Ntrials" %in% names(dat)) ".Ntrials" else NULL
    if (!is.null(nt_col)) inla_args$Ntrials <- dat[[nt_col]]
  }
  if (design$inla_family == "poisson" && ".E" %in% names(dat)) {
    inla_args$E <- dat$.E
  }

  fit <- tryCatch(
    suppressWarnings(suppressMessages(do.call(INLA::inla, inla_args))),
    error = function(e) e
  )
  if (inherits(fit, "error") || is.null(fit$summary.fixed)) {
    stop("INLA model fitting failed at this look: ",
         if (inherits(fit, "error")) conditionMessage(fit) else "no fixed effects returned.",
         call. = FALSE)
  }

  fitnames <- rownames(fit$summary.fixed)
  target <- if (design$effect_name %in% fitnames) {
    design$effect_name
  } else {
    cand <- grep(paste0("^", design$effect_name), fitnames, value = TRUE)
    if (length(cand) >= 1L) cand[1L] else NA_character_
  }
  if (is.na(target)) {
    stop("Monitored effect '", design$effect_name,
         "' not found among fitted coefficients: ",
         paste(fitnames, collapse = ", "), call. = FALSE)
  }

  mu  <- as.numeric(fit$summary.fixed[target, "mean"])
  sdv <- as.numeric(fit$summary.fixed[target, "sd"])
  z   <- stats::qnorm(1 - (1 - design$credible_level) / 2)
  ci_lower <- mu - z * sdv
  ci_upper <- mu + z * sdv

  pr_in <- NA_real_
  if (identical(design$metric, "direction")) {
    pr <- if (identical(design$alternative, "greater"))
      1 - stats::pnorm(0, mu, sdv) else stats::pnorm(0, mu, sdv)
  } else if (identical(design$metric, "threshold")) {
    pr <- if (identical(design$alternative, "greater"))
      1 - stats::pnorm(design$effect_threshold, mu, sdv)
    else stats::pnorm(design$effect_threshold, mu, sdv)
  } else {
    pr_in <- stats::pnorm(design$rope_bounds[2L], mu, sdv) -
      stats::pnorm(design$rope_bounds[1L], mu, sdv)
    pr <- 1 - pr_in
  }

  # ---- Decision --------------------------------------------------------------
  decision <- "continue"
  if (identical(design$metric, "rope")) {
    if (pr >= ps)            decision <- "stop_success"
    else if (pr_in >= ps)    decision <- "stop_equivalence"
  } else {
    if (pr >= ps)                          decision <- "stop_success"
    else if (!is.na(pf) && pr <= pf)       decision <- "stop_futility"
  }
  if (identical(decision, "continue") && n >= design$n_max) {
    decision <- "max_n_inconclusive"
  }

  row <- tibble::tibble(
    look        = k,
    time        = format(Sys.time(), tz = "UTC", usetz = TRUE),
    n           = n,
    planned_n   = planned_n,
    est_mean    = mu,
    est_sd      = sdv,
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    pr          = pr,
    pr_in_rope  = pr_in,
    threshold_success  = ps,
    threshold_futility = pf,
    decision    = decision,
    deviations  = if (length(deviations)) paste(deviations, collapse = "; ") else "",
    note        = note %||% ""
  )

  looks_df <- if (is.null(monitor)) row else dplyr::bind_rows(monitor$looks, row)
  status <- if (identical(decision, "continue")) "ongoing" else decision

  out <- list(
    design      = design,
    fingerprint = design$fingerprint,
    looks       = looks_df,
    status      = status
  )
  class(out) <- c("powerbrmsINLA_seq_monitor", "list")
  out
}


#' Print method for sequential analysis monitor objects
#'
#' @param x A `"powerbrmsINLA_seq_monitor"` object.
#' @param digits Digits for printed numbers.
#' @param ... Unused; present for S3 compatibility.
#' @return `x`, invisibly.
#' @export
print.powerbrmsINLA_seq_monitor <- function(x, digits = 3L, ...) {
  d <- x$design
  cat("Sequential Bayesian analysis monitor (powerbrmsINLA)\n")
  cat("=====================================================\n")
  cat("Study       :", d$label, "\n")
  cat("Design      :", d$metric, "metric; fingerprint", x$fingerprint, "\n")
  cat("Status      :", x$status, "\n\n")

  show <- as.data.frame(x$looks[, c("look", "n", "planned_n", "est_mean",
                                    "ci_lower", "ci_upper", "pr",
                                    "threshold_success", "decision")])
  num <- vapply(show, is.numeric, logical(1L))
  show[num] <- lapply(show[num], round, digits = digits)
  print(show, row.names = FALSE)

  dev_rows <- x$looks$deviations[nzchar(x$looks$deviations)]
  if (length(dev_rows) > 0L) {
    cat("\nRecorded deviations:\n")
    for (i in which(nzchar(x$looks$deviations))) {
      cat("  look", x$looks$look[i], ":", x$looks$deviations[i], "\n")
    }
  }

  last <- x$looks[nrow(x$looks), ]
  if (identical(x$status, "stop_success") && last$n < d$n_max) {
    cat("\nNote: the study stopped early for success. Effect estimates at\n")
    cat("early stops tend to be exaggerated; report the full posterior and\n")
    cat("the design's simulated operating characteristics (see\n")
    cat("brms_inla_sequential_trial()) alongside the decision.\n")
  }
  if (identical(x$status, "ongoing")) {
    nxt <- d$looks[min(nrow(x$looks) + 1L, d$n_looks)]
    cat("\nNext planned look at n =", nxt, "(of n_max =", paste0(d$n_max, ")."), "\n")
  }
  invisible(x)
}


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

#' Plot a sequential analysis trajectory with stopping boundaries
#'
#' @description
#' Visualises a [sequential_analysis()] monitor: either the monitored
#' posterior probability across looks with the prespecified success and
#' futility boundaries, or the effect estimate with its credible interval
#' across looks.
#'
#' @param monitor A `"powerbrmsINLA_seq_monitor"` object.
#' @param type `"probability"` (default) or `"estimate"`.
#' @param title,subtitle Optional plot annotations.
#' @return A ggplot object.
#' @export
plot_sequential_monitor <- function(
    monitor,
    type = c("probability", "estimate"),
    title = NULL,
    subtitle = NULL
) {
  if (!inherits(monitor, "powerbrmsINLA_seq_monitor")) {
    stop("`monitor` must be a 'powerbrmsINLA_seq_monitor' object from ",
         "sequential_analysis().", call. = FALSE)
  }
  type <- match.arg(type)
  d    <- monitor$design
  df   <- as.data.frame(monitor$looks)

  if (identical(type, "probability")) {
    bounds <- data.frame(
      n       = d$looks,
      success = d$prob_success,
      futility = if (is.null(d$prob_futility)) NA_real_ else d$prob_futility
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$n, y = .data$pr)) +
      ggplot2::geom_step(
        data = bounds,
        ggplot2::aes(x = .data$n, y = .data$success),
        inherit.aes = FALSE, colour = "#22A884", linetype = "dashed"
      )
    if (!is.null(d$prob_futility)) {
      p <- p + ggplot2::geom_step(
        data = bounds,
        ggplot2::aes(x = .data$n, y = .data$futility),
        inherit.aes = FALSE, colour = "#B85042", linetype = "dashed"
      )
    }
    p <- p +
      .geom_line_lw(width = 0.8) +
      .geom_point_lw(width = 2.5) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        title = title %||% paste0(d$label, ": monitored posterior probability"),
        subtitle = subtitle %||% paste0(
          "Metric: ", d$metric,
          "; dashed lines: success",
          if (!is.null(d$prob_futility)) " / futility" else "",
          " boundaries; fingerprint ", substr(monitor$fingerprint, 1, 8)
        ),
        x = "Accumulated sample size (n)",
        y = "Monitored posterior probability"
      ) +
      ggplot2::theme_minimal()
    return(p)
  }

  # type == "estimate"
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$n, y = .data$est_mean)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      fill = "#2A788E", alpha = 0.18
    ) +
    .geom_line_lw(width = 0.8) +
    .geom_point_lw(width = 2.5)
  if (identical(d$metric, "threshold")) {
    p <- p + ggplot2::geom_hline(yintercept = d$effect_threshold,
                                 linetype = "dotted")
  } else if (identical(d$metric, "rope")) {
    p <- p + ggplot2::geom_hline(yintercept = d$rope_bounds,
                                 linetype = "dotted")
  } else {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dotted")
  }
  p +
    ggplot2::labs(
      title = title %||% paste0(d$label, ": effect estimate across looks"),
      subtitle = subtitle %||% paste0(
        100 * d$credible_level, "% credible interval; fingerprint ",
        substr(monitor$fingerprint, 1, 8)
      ),
      x = "Accumulated sample size (n)",
      y = paste0("Posterior mean of '", d$effect_name, "'")
    ) +
    ggplot2::theme_minimal()
}

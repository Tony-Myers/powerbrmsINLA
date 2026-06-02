
#' Map a brms Family to an INLA Family (Modern, Robust)
#'
#' @param family brms family object or string (e.g., gaussian(), binomial())
#' @return List with $inla and $brms family names.
#' @keywords internal
.to_inla_family <- function(family) {
  # Step 1: Extract family name safely
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

  # Step 2: Map to INLA family name
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

  # Step 3: Return result
  return(list(inla = inla_family, brms = family_name))
}

#' Parse brms-like Random Effects Terms (Modern Robust)
#'
#' Extracts random-effects specifications from a brms-style formula string.
#' Supports (1 | g), (1 + x | g) - only one slope per term (warns otherwise!).
#' @param formula A model formula.
#' @return List of RE spec lists (group, intercept/slope info)
#' @keywords internal
.parse_re_terms <- function(formula) {
  ftxt <- paste(deparse(formula), collapse = "")
  raw <- gregexpr("\\([^\\|\\)]+\\|[^\\)]+\\)", ftxt, perl = TRUE)
  matches <- regmatches(ftxt, raw)[[1]]
  if (length(matches) == 0L) return(list())
  re_specs <- list(); idx <- 0L
  for (m in matches) {
    inside <- sub("^\\(", "", sub("\\)$", "", m))
    parts  <- strsplit(inside, "\\|", fixed = FALSE)[[1]]
    lhs    <- trimws(parts[1])
    grp    <- trimws(parts[2])
    comps  <- trimws(strsplit(lhs, "\\+")[[1]])
    has_intercept <- any(comps == "1")
    slope_vars    <- setdiff(comps, "1")
    if (length(slope_vars) > 1L) {
      warning("Only one random slope per RE term is supported. Ignoring extras: ",
              paste(slope_vars[-1], collapse = ", "))
      slope_vars <- slope_vars[1]
    }
    idx <- idx + 1L
    re_specs[[idx]] <- list(
      group        = grp,
      has_intercept= has_intercept,
      slope        = if (length(slope_vars) == 1L) slope_vars[[1]] else NULL,
      id_intercept = paste0(".re_id_", idx, "_int"),
      id_slope     = if (length(slope_vars) == 1L) paste0(".re_id_", idx, "_slope") else NULL
    )
  }
  re_specs
}

#' Convert brms Formula to INLA Formula (Multi-Fixed Support)
#'
#' Converts brms-style formula to INLA-ready formula, robust to multi-effect, and random effects.
#' @param formula brms-style formula.
#' @param drop_fixed Character vector of fixed effects to drop (optional).
#' @param hyper_by_re Optional named list of INLA hyper prior specifications by
#'   random-effect grouping factor.
#' @return List: $inla_formula, $re_specs, and $re_hyper_groups.
#' @keywords internal
.brms_to_inla_formula2 <- function(formula, drop_fixed = NULL, hyper_by_re = NULL) {
  re_specs <- .parse_re_terms(formula)
  re_hyper_groups <- character(0L)
  format_re_hyper <- function(hyper) {
    if (is.null(hyper) || is.null(hyper$prec)) return(NULL)
    prec <- hyper$prec
    if (is.null(prec$prior) || is.null(prec$param)) return(NULL)
    param <- paste(format(as.numeric(prec$param), digits = 15, scientific = FALSE),
                   collapse = ", ")
    sprintf("list(prec = list(prior = '%s', param = c(%s)))", prec$prior, param)
  }
  # Remove RE from fixed part
  ftxt <- paste(deparse(formula), collapse = "")
  ftxt_fixed <- gsub("\\([^\\|\\)]+\\|[^\\)]+\\)", "", ftxt, perl = TRUE)
  ftxt_fixed <- gsub("\\+\\s*\\+", "+", ftxt_fixed)
  ftxt_fixed <- gsub("~\\s*\\+", "~", ftxt_fixed)
  ftxt_fixed <- gsub("\\+\\s*$", "", ftxt_fixed)
  f_fixed <- as.formula(ftxt_fixed)
  resp <- as.character(f_fixed[[2L]])
  tt <- terms(f_fixed)
  rhs_terms <- attr(tt, "term.labels")
  has_intercept <- attr(tt, "intercept") == 1
  if (!is.null(drop_fixed)) {
    rhs_terms <- setdiff(rhs_terms, drop_fixed)
  }
  rhs <- character(0)
  rhs <- c(rhs, if (has_intercept) "1" else "-1")
  if (length(rhs_terms) > 0L) rhs <- c(rhs, rhs_terms)
  # Add random effects
  for (re in re_specs) {
    hyper_expr <- format_re_hyper(hyper_by_re[[re$group]])
    hyper_arg <- if (!is.null(hyper_expr)) paste0(", hyper = ", hyper_expr) else ""
    if (!is.null(hyper_expr)) re_hyper_groups <- union(re_hyper_groups, re$group)
    if (isTRUE(re$has_intercept))
      rhs <- c(rhs, sprintf("f(%s, model='iid'%s)", re$id_intercept, hyper_arg))
    if (!is.null(re$slope))
      rhs <- c(rhs, sprintf("f(%s, %s, model='iid'%s)", re$id_slope, re$slope, hyper_arg))
  }
  rhs_str <- if (length(rhs) > 0L) paste(rhs, collapse = " + ") else "1"
  inla_f <- as.formula(paste(resp, "~", rhs_str))
  list(inla_formula = inla_f, re_specs = re_specs, re_hyper_groups = re_hyper_groups)
}

#' Map brms Priors to INLA Priors (Multi-Fixed)
#'
#' Parses a brms prior specification, mapping supported priors to INLA controls.
#' Supports normal and student_t priors for fixed effects, and exponential
#' priors on residual/group-level SDs via INLA PC precision priors. Unsupported
#' priors are recorded in the audit instead of being silently translated.
#' @param priors brms prior object or data.frame.
#' @param approx_t_as_normal Logical; student_t priors with df > 2 are treated as normal.
#' @param pc_alpha Tail probability used for PC precision priors.
#' @param family_control_supplied Logical; if TRUE, translated sigma priors are
#'   audited but not applied because direct INLA family controls take precedence.
#' @param inla_family INLA family name, when known.
#' @return List with $control_fixed, $control_family, $hyper_by_re, and $prior_audit.
#' @keywords internal
.map_brms_priors_to_inla <- function(priors, approx_t_as_normal = TRUE,
                                     pc_alpha = 0.05,
                                     family_control_supplied = FALSE,
                                     inla_family = NULL) {
  out <- list(control_fixed = list(), control_family = list(), hyper_by_re = list(),
              prior_audit = .empty_prior_audit())
  if (is.null(priors)) return(out)
  pr_df <- tryCatch(as.data.frame(priors), error = function(e) NULL)
  if (is.null(pr_df) || !"class" %in% names(pr_df) || !"prior" %in% names(pr_df)) return(out)
  value_or_empty <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  }
  parse_normal <- function(s) {
    m <- regmatches(s, regexec("^\\s*normal\\s*\\(([^,]+),\\s*([^\\)]+)\\)\\s*$", s, ignore.case = TRUE))[[1]]
    if (length(m) == 3) list(mu = as.numeric(m[2]), sigma = as.numeric(m[3])) else NULL
  }
  parse_student_t <- function(s) {
    m <- regmatches(s, regexec("^\\s*student_t\\s*\\(([^,]+),\\s*([^,]+),\\s*([^\\)]+)\\)\\s*$", s, ignore.case = TRUE))[[1]]
    if (length(m) == 4) list(df = as.numeric(m[2]), mu = as.numeric(m[3]), sigma = as.numeric(m[4])) else NULL
  }
  parse_exponential <- function(s) {
    m <- regmatches(s, regexec("^\\s*exponential\\s*\\(([^\\)]+)\\)\\s*$", s, ignore.case = TRUE))[[1]]
    if (length(m) == 2) list(rate = as.numeric(m[2])) else NULL
  }
  pc_prec_from_exp_sd <- function(rate) {
    u <- -log(pc_alpha) / rate
    list(prec = list(prior = "pc.prec", param = c(u, pc_alpha)))
  }
  param_text <- function(x) paste0("c(", paste(format(as.numeric(x), digits = 15), collapse = ", "), ")")
  add_audit <- function(prior, cls, coef, group, target, inla_prior, inla_param,
                        method, applied, note) {
    out$prior_audit <<- rbind(
      out$prior_audit,
      data.frame(
        prior = prior, class = cls, coef = coef, group = group,
        target = target, inla_prior = inla_prior, inla_param = inla_param,
        method = method, applied = applied, note = note,
        stringsAsFactors = FALSE
      )
    )
  }
  mean_named <- list()
  prec_named <- list()
  mean_global <- NULL
  prec_global <- NULL
  mean_intercept <- NULL
  prec_intercept <- NULL
  cls_vec <- value_or_empty(pr_df$class)
  coef_vec <- if ("coef" %in% names(pr_df)) value_or_empty(pr_df$coef) else rep("", nrow(pr_df))
  group_vec <- if ("group" %in% names(pr_df)) value_or_empty(pr_df$group) else rep("", nrow(pr_df))
  for (i in seq_len(nrow(pr_df))) {
    cls <- cls_vec[i]
    coef <- coef_vec[i]
    group <- group_vec[i]
    pr_s <- value_or_empty(pr_df$prior[i])
    if (!is.character(pr_s) || !nzchar(pr_s)) next
    mu <- NA_real_; sigma <- NA_real_
    fixed_method <- NULL
    if (grepl("^\\s*normal\\s*\\(", pr_s, ignore.case = TRUE)) {
      p <- parse_normal(pr_s)
      if (!is.null(p)) { mu <- p$mu; sigma <- p$sigma; fixed_method <- "direct_normal_precision" }
    } else if (grepl("^\\s*student_t\\s*\\(", pr_s, ignore.case = TRUE)) {
      p <- parse_student_t(pr_s)
      if (!is.null(p)) {
        if (approx_t_as_normal && is.finite(p$df) && p$df > 2) {
          var_t <- (p$df / (p$df - 2)) * (p$sigma^2)
          mu <- p$mu; sigma <- sqrt(var_t); fixed_method <- "student_t_moment_matched_normal"
        } else {
          mu <- p$mu; sigma <- p$sigma; fixed_method <- "student_t_scale_as_normal"
        }
      }
    } else if (grepl("^\\s*exponential\\s*\\(", pr_s, ignore.case = TRUE)) {
      p <- parse_exponential(pr_s)
      if (is.null(p) || !is.finite(p$rate) || p$rate <= 0) {
        add_audit(pr_s, cls, coef, group, "", "", "", "unsupported", FALSE,
                  "Could not parse exponential(rate) with rate > 0.")
        next
      }
      hyper <- pc_prec_from_exp_sd(p$rate)
      if (identical(cls, "sigma")) {
        if (!is.null(inla_family) && !identical(inla_family, "gaussian")) {
          add_audit(pr_s, cls, coef, group, "control.family$hyper$prec",
                    "pc.prec", param_text(hyper$prec$param), "unsupported", FALSE,
                    "Residual SD prior translation is currently limited to Gaussian INLA models.")
        } else if (isTRUE(family_control_supplied)) {
          add_audit(pr_s, cls, coef, group, "control.family$hyper$prec",
                    "pc.prec", param_text(hyper$prec$param), "pc_prec_from_exponential_sd", FALSE,
                    "Not applied because user-supplied family_control takes precedence.")
        } else {
          out$control_family <- list(hyper = hyper)
          add_audit(pr_s, cls, coef, group, "control.family$hyper$prec",
                    "pc.prec", param_text(hyper$prec$param), "pc_prec_from_exponential_sd", TRUE,
                    "Translated exponential prior on residual SD to INLA PC precision prior.")
        }
      } else if (identical(cls, "sd")) {
        if (!nzchar(group)) {
          add_audit(pr_s, cls, coef, group, "", "pc.prec",
                    param_text(hyper$prec$param), "unsupported", FALSE,
                    "Group-level SD priors require a brms group name.")
        } else if (nzchar(coef)) {
          add_audit(pr_s, cls, coef, group, "", "pc.prec",
                    param_text(hyper$prec$param), "unsupported", FALSE,
                    "Coefficient-specific group-level SD priors are not translated by the minimal iid bridge.")
        } else {
          out$hyper_by_re[[group]] <- hyper
          add_audit(pr_s, cls, coef, group, paste0("f(... group = ", group, ")$hyper$prec"),
                    "pc.prec", param_text(hyper$prec$param), "pc_prec_from_exponential_sd", TRUE,
                    "Translated exponential prior on group-level SD to INLA PC precision prior for matching iid random effects.")
        }
      } else {
        add_audit(pr_s, cls, coef, group, "", "", "", "unsupported", FALSE,
                  "Only exponential priors on class 'sigma' and class 'sd' are translated.")
      }
      next
    } else {
      add_audit(pr_s, cls, coef, group, "", "", "", "unsupported", FALSE,
                "This brms prior family is not translated by the minimal INLA bridge.")
      next
    }
    if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
      add_audit(pr_s, cls, coef, group, "", "", "", "unsupported", FALSE,
                "Could not parse finite location and positive scale for this fixed-effect prior.")
      next
    }
    if (identical(cls, "Intercept")) {
      mean_intercept <- mu
      prec_intercept <- 1 / (sigma^2)
      add_audit(pr_s, cls, coef, group, "control.fixed$mean.intercept/prec.intercept",
                "normal", paste0("mean=", mu, "; precision=", prec_intercept),
                fixed_method, TRUE,
                "Direct INLA intercept prior. Note: brms intercept priors may differ because brms can center population-level predictors.")
    } else if (identical(cls, "b")) {
      if (!is.na(coef) && nzchar(coef)) {
        mean_named[[coef]] <- mu
        prec_named[[coef]] <- 1 / (sigma^2)
        add_audit(pr_s, cls, coef, group, paste0("control.fixed$mean/prec$", coef),
                  "normal", paste0("mean=", mu, "; precision=", 1 / (sigma^2)),
                  fixed_method, TRUE, "Direct fixed-effect Normal prior on the INLA linear predictor scale.")
      } else {
        mean_global <- mu
        prec_global <- 1 / (sigma^2)
        add_audit(pr_s, cls, coef, group, "control.fixed$mean/prec",
                  "normal", paste0("mean=", mu, "; precision=", 1 / (sigma^2)),
                  fixed_method, TRUE, "Global fixed-effect prior applied to non-intercept fixed effects.")
      }
    } else {
      add_audit(pr_s, cls, coef, group, "", "", "", "unsupported", FALSE,
                "Normal/student_t translation is currently limited to class 'b' and class 'Intercept'.")
    }
  }
  cf <- list()
  if (!is.null(mean_intercept)) cf$mean.intercept <- mean_intercept
  if (!is.null(prec_intercept)) cf$prec.intercept <- prec_intercept
  if (length(mean_named) > 0 && !is.null(mean_global)) {
    global_rows <- which(out$prior_audit$class == "b" & out$prior_audit$coef == "" &
                           out$prior_audit$target == "control.fixed$mean/prec")
    if (length(global_rows)) {
      out$prior_audit$applied[global_rows] <- FALSE
      out$prior_audit$method[global_rows] <- "unsupported"
      out$prior_audit$note[global_rows] <- paste0(
        out$prior_audit$note[global_rows],
        " Not applied because coefficient-specific fixed-effect priors were also supplied."
      )
    }
  }
  if (length(mean_named) > 0) cf$mean <- mean_named else if (!is.null(mean_global)) cf$mean <- mean_global
  if (length(prec_named) > 0) cf$prec <- prec_named else if (!is.null(prec_global)) cf$prec <- prec_global
  out$control_fixed <- cf
  out
}

.empty_prior_audit <- function() {
  data.frame(
    prior = character(), class = character(), coef = character(), group = character(),
    target = character(), inla_prior = character(), inla_param = character(),
    method = character(), applied = logical(), note = character(),
    stringsAsFactors = FALSE
  )
}

.mark_unmatched_re_priors <- function(prior_map, applied_groups) {
  audit <- prior_map$prior_audit
  if (is.null(audit) || !nrow(audit)) return(prior_map)
  rows <- which(audit$class == "sd" & audit$applied & !(audit$group %in% applied_groups))
  if (length(rows)) {
    audit$applied[rows] <- FALSE
    audit$method[rows] <- "unsupported"
    audit$note[rows] <- paste0(
      audit$note[rows],
      " Not applied because no matching iid random-effect term was generated from the formula."
    )
  }
  prior_map$prior_audit <- audit
  prior_map
}

.audit_re_correlation_terms <- function(prior_map, re_specs) {
  if (!length(re_specs)) return(prior_map)
  audit <- prior_map$prior_audit %||% .empty_prior_audit()
  for (re in re_specs) {
    if (!isTRUE(re$has_intercept) || is.null(re$slope)) next
    rows <- which(audit$class == "cor" & audit$group == re$group)
    note <- paste0(
      "Formula implies a correlated random intercept/slope term for group '",
      re$group,
      "'. The current INLA bridge represents this as independent iid terms; ",
      "correlation priors are not translated."
    )
    if (length(rows)) {
      audit$note[rows] <- paste(audit$note[rows], note)
      audit$applied[rows] <- FALSE
      audit$method[rows] <- "unsupported"
    } else {
      audit <- rbind(
        audit,
        data.frame(
          prior = "", class = "cor", coef = "", group = re$group,
          target = paste0("random-effect correlation for group ", re$group),
          inla_prior = "", inla_param = "", method = "unsupported",
          applied = FALSE, note = note,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  prior_map$prior_audit <- audit
  prior_map
}

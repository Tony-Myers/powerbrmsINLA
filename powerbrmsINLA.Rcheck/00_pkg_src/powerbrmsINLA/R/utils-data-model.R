
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
#' @return List: $inla_formula, $re_specs
#' @keywords internal
.brms_to_inla_formula2 <- function(formula, drop_fixed = NULL) {
  re_specs <- .parse_re_terms(formula)
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
    if (isTRUE(re$has_intercept))
      rhs <- c(rhs, sprintf("f(%s, model='iid')", re$id_intercept))
    if (!is.null(re$slope))
      rhs <- c(rhs, sprintf("f(%s, %s, model='iid')", re$id_slope, re$slope))
  }
  rhs_str <- if (length(rhs) > 0L) paste(rhs, collapse = " + ") else "1"
  inla_f <- as.formula(paste(resp, "~", rhs_str))
  list(inla_formula = inla_f, re_specs = re_specs)
}

#' Map brms Priors to INLA Priors (Multi-Fixed)
#'
#' Parses a brms prior specification, mapping fixed effect priors to INLA control.fixed.
#' Supports normal and student_t (variance-matched if df > 2).
#' @param priors brms prior object or data.frame.
#' @param approx_t_as_normal Logical; student_t priors with df > 2 are treated as normal.
#' @return List with $control_fixed and $hyper_by_re.
#' @keywords internal
.map_brms_priors_to_inla <- function(priors, approx_t_as_normal = TRUE) {
  out <- list(control_fixed = list(), hyper_by_re = list())
  if (is.null(priors)) return(out)
  pr_df <- tryCatch(as.data.frame(priors), error = function(e) NULL)
  if (is.null(pr_df) || !"class" %in% names(pr_df) || !"prior" %in% names(pr_df)) return(out)
  parse_normal <- function(s) {
    m <- regmatches(s, regexec("^\\s*normal\\s*\\(([^,]+),\\s*([^\\)]+)\\)\\s*$", s, ignore.case = TRUE))[[1]]
    if (length(m) == 3) list(mu = as.numeric(m[2]), sigma = as.numeric(m[3])) else NULL
  }
  parse_student_t <- function(s) {
    m <- regmatches(s, regexec("^\\s*student_t\\s*\\(([^,]+),\\s*([^,]+),\\s*([^\\)]+)\\)\\s*$", s, ignore.case = TRUE))[[1]]
    if (length(m) == 4) list(df = as.numeric(m[2]), mu = as.numeric(m[3]), sigma = as.numeric(m[4])) else NULL
  }
  mean_named <- list()
  prec_named <- list()
  mean_intercept <- NULL
  prec_intercept <- NULL
  for (i in seq_len(nrow(pr_df))) {
    cls <- pr_df$class[i]
    coef <- pr_df$coef[i]
    pr_s <- pr_df$prior[i]
    if (!is.character(pr_s) || !nzchar(pr_s)) next
    mu <- NA_real_; sigma <- NA_real_
    if (grepl("^\\s*normal\\s*\\(", pr_s, ignore.case = TRUE)) {
      p <- parse_normal(pr_s)
      if (!is.null(p)) { mu <- p$mu; sigma <- p$sigma }
    } else if (grepl("^\\s*student_t\\s*\\(", pr_s, ignore.case = TRUE)) {
      p <- parse_student_t(pr_s)
      if (!is.null(p)) {
        if (approx_t_as_normal && is.finite(p$df) && p$df > 2) {
          var_t <- (p$df / (p$df - 2)) * (p$sigma^2)
          mu <- p$mu; sigma <- sqrt(var_t)
        } else {
          mu <- p$mu; sigma <- p$sigma
        }
      }
    } else { next }
    if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) next
    if (identical(cls, "Intercept")) {
      mean_intercept <- mu
      prec_intercept <- 1 / (sigma^2)
    } else if (identical(cls, "b")) {
      if (!is.na(coef) && nzchar(coef)) {
        mean_named[[coef]] <- mu
        prec_named[[coef]] <- 1 / (sigma^2)
      }
    }
  }
  cf <- list()
  if (!is.null(mean_intercept)) cf$mean.intercept <- mean_intercept
  if (!is.null(prec_intercept)) cf$prec.intercept <- prec_intercept
  if (length(mean_named) > 0) cf$mean <- mean_named
  if (length(prec_named) > 0) cf$prec <- prec_named
  out$control_fixed <- cf
  out
}

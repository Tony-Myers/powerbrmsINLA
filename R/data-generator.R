#' Automatic Data Generator for brms + INLA Simulation (Multi-Effect Ready)
#'
#' Creates a simulation function taking n (sample size) and a named effect vector/list.
#' Supports multiple predictors, brms-style random effects, and most GLM families.
#'
#' @param formula Model formula (may include brms-style random effects).
#' @param effect_name Character vector of fixed effect names to manipulate.
#' @param family brms family object, e.g. gaussian(), binomial().
#' @param family_args Named list of family-specific arguments.
#' @param error_sd Residual SD for Gaussian-like families.
#' @param group_sd SD of random effects.
#' @param obs_per_group Number of observations per grouping level.
#' @param predictor_means Named list of means for continuous predictors.
#' @param predictor_sds Named list of SDs for continuous predictors.
#' @return Function: function(n, effect) returns a data.frame with n rows.
#' @keywords internal
.auto_data_generator <- function(
    formula,
    effect_name,
    family = gaussian(),
    family_args = list(),
    error_sd = 1,
    group_sd = 0.5,
    obs_per_group = 10,
    predictor_means = NULL,
    predictor_sds = NULL
) {
  re_specs <- .parse_re_terms(formula)
  ftxt <- paste(deparse(formula), collapse = "")
  ftxt_fixed <- gsub("\\([^\\|\\)]+\\|[^\\)]+\\)", "", ftxt, perl = TRUE)
  ftxt_fixed <- gsub("\\+\\s*\\+", "+", ftxt_fixed)
  ftxt_fixed <- gsub("~\\s*\\+", "~", ftxt_fixed)
  ftxt_fixed <- gsub("\\+\\s*$", "", ftxt_fixed)
  f_fixed <- as.formula(ftxt_fixed)
  resp <- as.character(f_fixed[[2L]])
  tt <- terms(f_fixed)
  fixed_terms <- attr(tt, "term.labels")
  all_vars <- if (length(fixed_terms)) unique(unlist(strsplit(fixed_terms, ":", fixed = TRUE))) else character(0)
  predictor_vars <- setdiff(all_vars, resp)
  
  fam_map <- .to_inla_family(family)
  fam_inla <- fam_map$inla
  
  # Family-specific defaults
  nb_theta <- family_args$theta %||% 1.0
  beta_phi <- family_args$phi %||% 50
  bb_phi <- family_args$phi %||% 50
  bb_size <- family_args$size %||% 10L
  E_exposure <- family_args$E %||% NULL
  t_df <- family_args$df %||% 5
  sn_alpha <- family_args$alpha %||% 3
  vm_kappa <- family_args$kappa %||% 2
  gev_xi <- family_args$xi %||% 0.1
  gev_sigma <- family_args$sigma %||% 1.0
  
  clamp01 <- function(x, eps = 1e-6) pmin(1 - eps, pmax(eps, x))
  
  # --- The returned function ---
  function(n, effect) {
    # effect is either named vector/list (multi-predictor) or scalar (single-predictor)
    if (is.list(effect) && !is.null(names(effect))) effect <- unlist(effect, use.names=TRUE)
    
    n_groups_list <- list()
    if (length(re_specs) > 0L) {
      for (i in seq_along(re_specs)) {
        group_name <- re_specs[[i]]$group
        n_groups_list[[group_name]] <- max(3L, round(n / obs_per_group))
      }
    }
    dat <- data.frame(row.names = seq_len(n))
    
    # Add grouping factors
    if (length(n_groups_list) > 0L) {
      for (group_name in names(n_groups_list)) {
        dat[[group_name]] <- factor(rep(seq_len(n_groups_list[[group_name]]), length.out = n))
      }
    }
    
    # Add predictors
    for (var in predictor_vars) {
      if (grepl("treatment|group|condition", var, ignore.case = TRUE)) {
        dat[[var]] <- factor(sample(c(0, 1), n, replace = TRUE))
      } else if (grepl("time|wave|visit", var, ignore.case = TRUE)) {
        if (length(n_groups_list) > 0L) {
          subject_var <- names(n_groups_list)[1L]
          n_subjects <- n_groups_list[[subject_var]]
          times_per_subject <- max(1L, round(n / n_subjects))
          dat[[var]] <- rep(0:(times_per_subject - 1L), n_subjects)[1L:n]
        } else {
          dat[[var]] <- seq_len(n) - 1L
        }
      } else {
        var_mean <- if (!is.null(predictor_means) && !is.null(predictor_means[[var]])) predictor_means[[var]] else 0
        var_sd <- if (!is.null(predictor_sds) && !is.null(predictor_sds[[var]])) predictor_sds[[var]] else 1
        dat[[var]] <- rnorm(n, mean = var_mean, sd = var_sd)
      }
    }
    
    # Linear predictor
    eta <- 0
    for (term in fixed_terms) {
      term_value <- if (grepl(":", term, fixed = TRUE)) {
        vars <- strsplit(term, ":", fixed = TRUE)[[1]]
        tv <- dat[[vars[1]]]
        for (v in vars[-1]) tv <- tv * if (is.factor(dat[[v]])) (as.numeric(dat[[v]]) - 1) else dat[[v]]
        tv
      } else {
        if (is.factor(dat[[term]])) (as.numeric(dat[[term]]) - 1) else dat[[term]]
      }
      # If this term matches effect_name, set its coefficient by value from 'effect'
      if (term %in% effect_name) {
        eta <- eta + effect[[term]] * term_value
      } else {
        coef_val <- rnorm(1L, 0, 0.3)  # default "noise" coefficient
        eta <- eta + coef_val * term_value
      }
    }
    # Random effects
    if (length(re_specs) > 0L) {
      for (re in re_specs) {
        glev <- n_groups_list[[re$group]]
        ri <- rnorm(glev, 0, group_sd)
        eta <- eta + ri[as.integer(dat[[re$group]])]
        if (!is.null(re$slope)) {
          rs <- rnorm(glev, 0, group_sd * 0.5)
          sval <- if (is.factor(dat[[re$slope]])) (as.numeric(dat[[re$slope]]) - 1) else dat[[re$slope]]
          eta <- eta + rs[as.integer(dat[[re$group]])] * sval
        }
      }
    }
    
    # Response generation
    if (fam_inla == "gaussian") {
      y <- eta + rnorm(n, 0, error_sd)
    } else if (fam_inla == "T") {
      y <- eta + error_sd * stats::rt(n, df = t_df)
    } else if (fam_inla == "binomial") {
      p <- clamp01(stats::plogis(eta))
      size <- as.integer(family_args$size %||% 1L)
      dat$.Ntrials <- size
      y <- stats::rbinom(n, size = size, prob = p)
    } else if (fam_inla == "poisson") {
      lam <- exp(eta)
      if (!is.null(E_exposure)) {
        dat$.E <- as.numeric(E_exposure)
        if (length(dat$.E) == 1L) dat$.E <- rep(dat$.E, n)
        lam <- lam * dat$.E
      }
      y <- stats::rpois(n, lambda = lam)
    } else if (fam_inla == "nbinomial") {
      mu <- exp(eta)
      y <- MASS::rnegbin(n, mu = mu, theta = nb_theta)
    } else if (fam_inla == "beta") {
      m <- clamp01(stats::plogis(eta))
      phi <- beta_phi
      y <- stats::rbeta(n, shape1 = m * phi, shape2 = (1 - m) * phi)
      y <- clamp01(y)
    } else if (fam_inla == "betabinomial") {
      m <- clamp01(stats::plogis(eta))
      phi <- bb_phi
      a <- m * phi; b <- (1 - m) * phi
      p <- stats::rbeta(n, a, b)
      size <- as.integer(family_args$size %||% bb_size)
      dat$.Ntrials <- size
      y <- stats::rbinom(n, size = size, prob = p)
    } else if (fam_inla == "lognormal") {
      y <- exp(eta + stats::rnorm(n, 0, error_sd))
    } else if (fam_inla == "weibull") {
      k <- family_args$shape %||% 1.5
      y <- stats::rweibull(n, shape = k, scale = exp(eta))
    } else if (fam_inla == "exponential") {
      y <- stats::rexp(n, rate = exp(-eta))
    } else if (fam_inla == "sn") {
      if (requireNamespace("sn", quietly = TRUE)) {
        y <- sn::rsn(n, xi = eta, omega = error_sd, alpha = sn_alpha)
      } else {
        warning("Package 'sn' not available; simulating skew-normal as Gaussian.")
        y <- eta + stats::rnorm(n, 0, error_sd)
      }
    } else if (fam_inla == "vm") {
      if (requireNamespace("circular", quietly = TRUE)) {
        mu <- circular::circular(eta %% (2 * pi))
        y <- as.numeric(circular::rvonmises(n, mu = mu, kappa = vm_kappa))
      } else {
        warning("Package 'circular' not available; simulating von Mises as wrapped Normal.")
        y <- (eta + stats::rnorm(n, 0, 1 / sqrt(vm_kappa + 1e-8))) %% (2 * pi)
      }
    } else if (fam_inla == "gev") {
      u <- stats::runif(n)
      xi <- gev_xi; sigma <- gev_sigma; mu <- eta
      if (abs(xi) < 1e-8) {
        y <- mu - sigma * log(-log(u))
      } else {
        y <- mu + (sigma / xi) * ((-log(u))^(-xi) - 1)
      }
    } else {
      stop("Family '", fam_inla, "' is not implemented in the generator.")
    }
    
    dat[[resp]] <- y
    
    # Add INLA RE index columns
    if (length(re_specs) > 0L) {
      for (re in re_specs) {
        if (isTRUE(re$has_intercept) && !is.null(re$id_intercept)) {
          dat[[re$id_intercept]] <- as.integer(dat[[re$group]])
        }
        if (!is.null(re$slope) && !is.null(re$id_slope)) {
          dat[[re$id_slope]] <- as.integer(dat[[re$group]])
        }
      }
    }
    
    dat
  }
}

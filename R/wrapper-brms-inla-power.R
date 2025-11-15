#' Design-based wrapper for Bayesian power / assurance
#'
#' Dispatches to one of the three engines depending on `design`.
#' This function must accept `...` and pass it on unchanged.
#'
#' @param design Character scalar: "fixed", "two_stage", or "sequential".
#' @param ... Arguments passed on to the corresponding engine.
#'
#' @return Whatever the underlying engine returns.
#'
#' @export
brms_inla_power_design <- function(
    design = c("fixed", "two_stage", "sequential"),
    ...
) {
  design <- match.arg(design)
  
  if (identical(design, "fixed")) {
    # Fixed-n engine from engine-main.R
    return(brms_inla_power(...))
  }
  
  if (identical(design, "two_stage")) {
    # Two-stage engine from engine-two-stage.R
    return(brms_inla_power_two_stage(...))
  }
  
  if (identical(design, "sequential")) {
    # Sequential engine from engine-sequential.R
    return(brms_inla_power_sequential(...))
  }
  
  stop("Unknown design: ", design, call. = FALSE)
}

#' Parallel wrapper for fixed-design Bayesian power / assurance simulations
#'
#' Parallelises over cells defined by sample_sizes x effect_grid for the
#' fixed-n engine brms_inla_power().
#'
#' @param design Character scalar. Currently only "fixed" is supported.
#' @param sample_sizes Numeric vector of sample sizes (required).
#' @param effect_grid Numeric vector or data frame defining effect scenarios
#'   (required).
#' @param nsims Integer number of simulations per cell.
#' @param n_cores Integer number of worker processes. Default is
#'   max(1L, parallel::detectCores() - 1L).
#' @param seed Integer base seed. Each cell uses seed + cell_id.
#' @param progress Logical or character; controls wrapper-level progress bar.
#' @param ... Further arguments passed directly to brms_inla_power(), such as
#'   formula, family, priors, effect_name, compute_bayes_factor, bf_method,
#'   inla_hyper, inla_num_threads, etc.
#'
#' @return A list with components summary, results, and settings.
#'
#' @export
brms_inla_power_parallel <- function(
    design       = c("fixed"),
    sample_sizes,
    effect_grid,
    nsims,
    n_cores      = max(1L, parallel::detectCores() - 1L),
    seed         = 123L,
    progress     = c("auto", "text", "none"),
    ...
) {
  ## design is only used for a sanity check, never forwarded
  design <- match.arg(design)
  if (!identical(design, "fixed")) {
    stop(
      "brms_inla_power_parallel currently supports only design = 'fixed'. ",
      "Use the sequential or two-stage engines directly for other designs.",
      call. = FALSE
    )
  }
  
  if (missing(sample_sizes)) stop("Argument 'sample_sizes' must be supplied.", call. = FALSE)
  if (missing(effect_grid))  stop("Argument 'effect_grid' must be supplied.",  call. = FALSE)
  if (missing(nsims))        stop("Argument 'nsims' must be supplied.",        call. = FALSE)
  
  ## Capture all extra args (formula, family, priors, bf_method, etc.)
  extra_args <- list(...)
  ## Ensure any 'design' that arrived via ... is dropped and never reaches the engine
  if ("design" %in% names(extra_args)) {
    extra_args$design <- NULL
  }
  
  ## Normalise progress
  if (is.logical(progress)) {
    progress_arg <- if (isTRUE(progress)) "text" else "none"
  } else {
    progress_arg <- match.arg(progress)
  }
  
  ## Sanitize n_cores
  if (!is.numeric(n_cores) || length(n_cores) != 1L || !is.finite(n_cores)) {
    n_cores <- 1L
  } else {
    n_cores <- as.integer(n_cores)
    if (n_cores < 1L) n_cores <- 1L
  }
  
  ## Sequential path: just call brms_inla_power() once
  if (n_cores == 1L) {
    return(
      do.call(
        brms_inla_power,
        c(
          list(
            sample_sizes = sample_sizes,
            effect_grid  = effect_grid,
            nsims        = nsims,
            seed         = seed,
            progress     = progress_arg
          ),
          extra_args
        )
      )
    )
  }
  
  ## Parallel path: split into cells
  is_multi_effect <- is.data.frame(effect_grid)
  effect_indices  <- if (is_multi_effect) seq_len(nrow(effect_grid)) else seq_along(effect_grid)
  n_indices       <- seq_along(sample_sizes)
  
  cell_grid <- expand.grid(
    n_index      = n_indices,
    effect_index = effect_indices,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  n_cells <- nrow(cell_grid)
  
  ## Per-cell function: calls brms_inla_power() directly
  run_cell <- function(cell_id) {
    idx <- cell_grid[cell_id, , drop = FALSE]
    single_n <- sample_sizes[idx$n_index]
    
    single_effect <- if (is_multi_effect) {
      effect_grid[idx$effect_index, , drop = FALSE]
    } else {
      effect_grid[idx$effect_index]
    }
    
    cell_seed <- seed + as.integer(cell_id)
    
    do.call(
      brms_inla_power,
      c(
        list(
          sample_sizes = single_n,
          effect_grid  = single_effect,
          nsims        = nsims,
          seed         = cell_seed,
          progress     = "none"
        ),
        extra_args
      )
    )
  }
  
  ## Create cluster
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  ## Export objects needed inside run_cell
  parallel::clusterExport(
    cl,
    varlist = c(
      "cell_grid",
      "sample_sizes",
      "effect_grid",
      "is_multi_effect",
      "seed",
      "nsims",
      "extra_args",
      "brms_inla_power"
    ),
    envir = environment()
  )
  
  ## Progress bar in parallel if pbapply is available
  use_pb <- (progress_arg %in% c("text", "auto")) &&
    requireNamespace("pbapply", quietly = TRUE)
  
  if (use_pb) {
    res_list <- pbapply::pblapply(
      X   = seq_len(n_cells),
      FUN = run_cell,
      cl  = cl
    )
  } else {
    res_list <- parallel::parLapply(cl, seq_len(n_cells), run_cell)
  }
  
  ## Combine summaries and results
  summaries <- lapply(res_list, function(x) x$summary)
  results   <- lapply(res_list, function(x) x$results)
  
  all_summary <- do.call(rbind, summaries)
  all_results <- do.call(rbind, results)
  
  rownames(all_summary) <- NULL
  rownames(all_results) <- NULL
  
  ## Use the first result as template
  out <- res_list[[1L]]
  
  out$summary <- all_summary
  out$results <- all_results
  
  if (is.null(out$settings)) {
    out$settings <- list()
  }
  if (!is.list(out$settings)) {
    out$settings <- list(settings = out$settings)
  }
  
  out$settings$design       <- "fixed"
  out$settings$effect_grid  <- effect_grid
  out$settings$sample_sizes <- sample_sizes
  out$settings$nsims        <- nsims
  out$settings$parallel     <- list(
    enabled  = TRUE,
    n_cores  = n_cores,
    n_cells  = n_cells,
    progress = if (use_pb) progress_arg else "none"
  )
  
  out
}
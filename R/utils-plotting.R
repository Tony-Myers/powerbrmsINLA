# Internal: get a results data frame from either a results list or a data frame
.get_results_df <- function(x) {
  if (is.list(x) && !is.null(x$results)) {
    return(x$results)
  }
  if (is.data.frame(x)) {
    return(x)
  }
  stop("Expected a list with $results or a data.frame.", call. = FALSE)
}

# Internal: get a summary data frame from either a results list or a data frame
.get_summary_df <- function(x) {
  if (is.list(x) && !is.null(x$summary)) {
    return(x$summary)
  }
  if (is.data.frame(x)) {
    return(x)
  }
  stop("Expected a list with $summary or a data.frame.", call. = FALSE)
}

# Internal: apply exact-match filters like list(treatment = 0.3)
.apply_effect_filters <- function(df, filters = NULL) {
  if (is.null(filters) || length(filters) == 0L) return(df)
  for (nm in names(filters)) {
    if (!nm %in% names(df)) {
      warning("Filter '", nm, "' not found in data; ignoring.", call. = FALSE)
      next
    }
    df <- df[df[[nm]] == filters[[nm]], , drop = FALSE]
  }
  df
}
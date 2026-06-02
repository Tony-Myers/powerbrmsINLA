pkgname <- "powerbrmsINLA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('powerbrmsINLA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("assurance_prior_weights")
### * assurance_prior_weights

flush(stderr()); flush(stdout())

### Name: assurance_prior_weights
### Title: Create prior weights over an effect grid for use with
###   compute_assurance()
### Aliases: assurance_prior_weights

### ** Examples

effects <- c(0.1, 0.3, 0.5, 0.7, 0.9)
# Normal prior centred on 0.5
assurance_prior_weights(effects, dist = "normal", mean = 0.5, sd = 0.2)
# Uniform prior (equal weight)
assurance_prior_weights(effects, dist = "uniform")



cleanEx()
nameEx("brms_inla_power")
### * brms_inla_power

flush(stderr()); flush(stdout())

### Name: brms_inla_power
### Title: Core Bayesian Assurance / Power Simulation (Modern, Multi-Effect
###   Ready)
### Aliases: brms_inla_power

### ** Examples

## Not run: 
##D # Integrate over uncertainty in the residual SD using a half-normal prior
##D # centred at 1.0 with spread 0.3
##D brms_inla_power(
##D   formula      = y ~ treatment,
##D   effect_name  = "treatment",
##D   effect_grid  = 0.5,
##D   sample_sizes = c(50, 100),
##D   nsims        = 50,
##D   error_sd     = list(dist = "halfnormal", sd = 0.3, location = 1.0),
##D   seed         = 42
##D )
## End(Not run)




cleanEx()
nameEx("brms_inla_power_sequential")
### * brms_inla_power_sequential

flush(stderr()); flush(stdout())

### Name: brms_inla_power_sequential
### Title: Sequential Bayesian Assurance Simulation Engine (Modern,
###   Multi-Effect Ready)
### Aliases: brms_inla_power_sequential

### ** Examples

## Not run: 
##D # Sequential design with automatic threading
##D results <- brms_inla_power_sequential(
##D   formula = outcome ~ treatment,
##D   effect_name = "treatment",
##D   effect_grid = c(0.2, 0.5, 0.8),
##D   sample_sizes = c(50, 100, 200),
##D   metric = "direction",
##D   target = 0.80
##D )
##D print(results$summary)
## End(Not run)



cleanEx()
nameEx("brms_inla_power_two_stage")
### * brms_inla_power_two_stage

flush(stderr()); flush(stdout())

### Name: brms_inla_power_two_stage
### Title: Two-Stage Bayesian Assurance Simulation (Multi-Effect,
###   User-Friendly API)
### Aliases: brms_inla_power_two_stage

### ** Examples

## Not run: 
##D # Two-stage design with threading
##D effect_grid <- expand.grid(
##D   treatment = c(0.2, 0.5, 0.8),
##D   covariate = c(0.1, 0.3)
##D )
##D results <- brms_inla_power_two_stage(
##D   formula = outcome ~ treatment + covariate,
##D   effect_name = c("treatment", "covariate"),
##D   effect_grid = effect_grid,
##D   n_range = c(50, 200),
##D   stage1_nsims = 3,
##D   stage2_nsims = 3,
##D    error_sd = 1 
##D )
##D print(results$summary)
## End(Not run)



cleanEx()
nameEx("compute_assurance")
### * compute_assurance

flush(stderr()); flush(stdout())

### Name: compute_assurance
### Title: Compute unconditional Bayesian assurance from simulation results
### Aliases: compute_assurance

### ** Examples

# Build a small synthetic power_result without running INLA
syn_summary <- data.frame(
  n               = rep(c(50, 100, 200), each = 3),
  treatment       = rep(c(0.2, 0.5, 0.8), 3),
  power_direction = c(0.40, 0.65, 0.85,
                      0.60, 0.82, 0.95,
                      0.72, 0.90, 0.98),
  stringsAsFactors = FALSE
)
syn_result <- list(
  summary  = syn_summary,
  settings = list(effect_name = "treatment")
)

# (a) Uniform weights — assurance is the simple mean of per-cell powers
w_uniform <- c("0.2" = 1/3, "0.5" = 1/3, "0.8" = 1/3)
out <- compute_assurance(syn_result, prior_weights = w_uniform)
print(out)

# (b) Normal design prior centred on a medium-sized effect
out2 <- compute_assurance(
  syn_result,
  prior_weights = list(dist = "normal", mean = 0.5, sd = 0.2)
)
print(out2)

# (c) Using assurance_prior_weights() to build the weight vector explicitly
w_norm <- assurance_prior_weights(c(0.2, 0.5, 0.8), dist = "normal",
                                  mean = 0.5, sd = 0.2)
out3 <- compute_assurance(syn_result, prior_weights = w_norm)



cleanEx()
nameEx("decide_sample_size")
### * decide_sample_size

flush(stderr()); flush(stdout())

### Name: decide_sample_size
### Title: Decide recommended sample size from power/assurance results
### Aliases: decide_sample_size

### ** Examples

# Build a small synthetic power_result without running INLA
syn_summary <- data.frame(
  n               = rep(c(50, 100, 200), each = 3),
  treatment       = rep(c(0.2, 0.5, 0.8), 3),
  power_direction = c(0.40, 0.65, 0.85,
                      0.60, 0.82, 0.95,
                      0.72, 0.90, 0.98),
  power_threshold = c(0.30, 0.55, 0.75,
                      0.50, 0.72, 0.88,
                      0.62, 0.80, 0.92),
  stringsAsFactors = FALSE
)
syn_result <- list(
  summary  = syn_summary,
  settings = list(effect_name = "treatment")
)

# --- Assurance mode: each metric value IS the assurance target ---
w <- assurance_prior_weights(c(0.2, 0.5, 0.8), dist = "normal",
                              mean = 0.5, sd = 0.2)
# Find n where direction assurance >= 0.80 AND threshold assurance >= 0.75
rec_assurance <- decide_sample_size(
  syn_result,
  direction     = 0.80,
  threshold     = 0.75,
  prior_weights = w
)
print(rec_assurance)

# --- Conditional mode (backward compatible) ---
rec_cond <- decide_sample_size(syn_result, direction = 0.80)
print(rec_cond)



cleanEx()
nameEx("plot_assurance_curve")
### * plot_assurance_curve

flush(stderr()); flush(stdout())

### Name: plot_assurance_curve
### Title: Plot Assurance Curve(s) vs Sample Size
### Aliases: plot_assurance_curve

### ** Examples

syn_summary <- data.frame(
  n               = rep(c(50, 100, 200), each = 3),
  treatment       = rep(c(0.2, 0.5, 0.8), 3),
  power_direction = c(0.40, 0.65, 0.85, 0.60, 0.82, 0.95, 0.72, 0.90, 0.98)
)
pr <- list(summary = syn_summary, settings = list(effect_name = "treatment"))
a1 <- compute_assurance(pr, list(dist = "normal", mean = 0.5, sd = 0.15))
a2 <- compute_assurance(pr, list(dist = "normal", mean = 0.5, sd = 0.30))
plot_assurance_curve(list(a1, a2), labels = c("Informative", "Diffuse"))




cleanEx()
nameEx("plot_design_prior")
### * plot_design_prior

flush(stderr()); flush(stdout())

### Name: plot_design_prior
### Title: Plot Design Prior Density over the Effect Grid
### Aliases: plot_design_prior

### ** Examples

plot_design_prior(
  prior_spec  = list(dist = "normal", mean = 0.5, sd = 0.15),
  effect_grid = c(0.2, 0.5, 0.8)
)

# Multiple priors
plot_design_prior(
  prior_spec  = list(
    list(dist = "normal", mean = 0.5, sd = 0.10),
    list(dist = "normal", mean = 0.5, sd = 0.30)
  ),
  effect_grid = c(0.2, 0.5, 0.8),
  labels      = c("Informative", "Diffuse")
)




cleanEx()
nameEx("plot_power_assurance_overlay")
### * plot_power_assurance_overlay

flush(stderr()); flush(stdout())

### Name: plot_power_assurance_overlay
### Title: Plot Conditional Power Curves with Assurance Overlay
### Aliases: plot_power_assurance_overlay

### ** Examples

syn_summary <- data.frame(
  n               = rep(c(50, 100, 200), each = 3),
  treatment       = rep(c(0.2, 0.5, 0.8), 3),
  power_direction = c(0.40, 0.65, 0.85, 0.60, 0.82, 0.95, 0.72, 0.90, 0.98)
)
pr <- list(summary = syn_summary, settings = list(effect_name = "treatment"))
a  <- compute_assurance(pr, list(dist = "normal", mean = 0.5, sd = 0.15))
plot_power_assurance_overlay(pr, a)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

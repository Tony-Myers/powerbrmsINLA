pkgname <- "powerbrmsINLA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "powerbrmsINLA-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('powerbrmsINLA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("brms_inla_power")
### * brms_inla_power

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: brms_inla_power
### Title: Core Bayesian Assurance / Power Simulation (Modern, Multi-Effect
###   Ready)
### Aliases: brms_inla_power

### ** Examples

## No test: 
# Basic usage with automatic INLA threading
results <- brms_inla_power(
  formula = outcome ~ treatment,
  effect_name = "treatment",
  effect_grid = c(0.2, 0.5, 0.8),
  sample_sizes = c(50, 100, 200),
  nsims = 3
)
print(results$summary)

# Manual INLA threading control
results <- brms_inla_power(
  formula = outcome ~ treatment,
  effect_name = "treatment",
  effect_grid = c(0.2, 0.5, 0.8),
  sample_sizes = c(50, 100, 200),
  inla_num_threads = "8:1",  # Use 8 threads for faster computation
  nsims = 3
)

# Multi-effect design with threading
effect_grid <- expand.grid(
  treatment = c(0, 0.3, 0.6),
  age_effect = c(0, 0.2)
)
results <- brms_inla_power(
  formula = outcome ~ treatment + age_effect,
  effect_name = c("treatment", "age_effect"),
  effect_grid = effect_grid,
  sample_sizes = c(100, 200, 400),
  nsims = 3
)
print(results$summary)
## End(No test)
# Quick parameter check (runs instantly)
formals(brms_inla_power)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("brms_inla_power", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("brms_inla_power_sequential")
### * brms_inla_power_sequential

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: brms_inla_power_sequential
### Title: Sequential Bayesian Assurance Simulation Engine (Modern,
###   Multi-Effect Ready)
### Aliases: brms_inla_power_sequential

### ** Examples

## No test: 
# Sequential design with automatic threading
results <- brms_inla_power_sequential(
  formula = outcome ~ treatment,
  effect_name = "treatment",
  effect_grid = c(0.2, 0.5, 0.8),
  sample_sizes = c(50, 100, 200),
  metric = "direction",
  target = 0.80
)
print(results$summary)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("brms_inla_power_sequential", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("brms_inla_power_two_stage")
### * brms_inla_power_two_stage

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: brms_inla_power_two_stage
### Title: Two-Stage Bayesian Assurance Simulation (Multi-Effect,
###   User-Friendly API)
### Aliases: brms_inla_power_two_stage

### ** Examples

## No test: 
# Two-stage design with threading
effect_grid <- expand.grid(
  treatment = c(0.2, 0.5, 0.8),
  covariate = c(0.1, 0.3)
)
results <- brms_inla_power_two_stage(
  formula = outcome ~ treatment + covariate,
  effect_name = c("treatment", "covariate"),
  effect_grid = effect_grid,
  n_range = c(50, 200),
  stage1_nsims = 3,
  stage2_nsims = 3,
   error_sd = 1 
)
print(results$summary)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("brms_inla_power_two_stage", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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

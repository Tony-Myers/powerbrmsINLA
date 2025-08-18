
# Test setup file for brmsINLApower package

# Suppress warnings during testing
options(warn = -1)

# Set up testing environment
library(testthat)
library(brms)
library(ggplot2)
library(dplyr)
library(tibble)

# Helper function for tests - FIX: Define expect_function
expect_ggplot <- function(object) {
  expect_s3_class(object, "ggplot")
}

expect_function <- function(object) {
  expect_type(object, "closure")
}

# Create standard test data for reuse
create_test_data <- function(n = 100) {
  data.frame(
    y = rnorm(n),
    treatment = factor(sample(c(0, 1), n, replace = TRUE)),
    age = rnorm(n),
    subject = factor(rep(1:(n/5), each = 5)[1:n]),
    time = rep(0:4, length.out = n)
  )
}

message("Test setup complete")

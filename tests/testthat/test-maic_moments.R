# maic_moments

library(testthat)
library(dplyr)

# Helper Generate Mock Data for Tests
generate_mock_data <- function() {
  set.seed(123)
  
  # Mock IPD
  ipd <- data.frame(
    trt = as.factor(sample(c("A", "C"), 100, replace = TRUE)),
    AGE = rnorm(100, mean = 50, sd = 8),
    SCORE = rnorm(100, mean = 10, sd = 2),
    y = rbinom(100, 1, 0.5)
  )
  
  # Mock ALD (Means and SDs slightly shifted from IPD)
  ald <- data.frame(
    variable = c("AGE", "AGE", "SCORE", "SCORE", "AGE:SCORE"),
    statistic = c("mean", "sd", "mean", "sd", "mean"),
    value = c(52.0, 9.5, 11.0, 2.5, 572.0) # 52 * 11 = 572 for the interaction
  )
  
  list(ipd = ipd, ald = ald)
}

# Tests for Moments = 2 (Variance Matching)
test_that("maic.boot successfully runs with moments = 2 when SD is provided", {
  mock <- generate_mock_data()
  
  # Run the bootstrap function manually for 1 iteration
  result <- outstandR:::maic.boot(
    ipd = mock$ipd,
    indices = 1:nrow(mock$ipd),
    outcome_model = y ~ trt,
    balance_model = ~ AGE + SCORE,
    family = binomial(),
    ald = mock$ald,
    trt_var = "trt",
    moments = 2,
    int = FALSE
  )
  
  # Check that it returns the expected 4-part vector structure
  expect_type(result, "double")
  expect_true("pC" %in% names(result))
  expect_true("pA" %in% names(result))
  expect_true("ESS" %in% names(result))
  
  # ESS should be a valid number greater than 1, but less than or equal to N
  ess_val <- result["ESS"]
  expect_true(is.numeric(ess_val))
  expect_true(ess_val > 1 && ess_val <= 100)
})

test_that("maic.boot throws an error if moments = 2 but SD is missing in ALD", {
  mock <- generate_mock_data()
  
  # Remove the SD for AGE to trigger the error
  bad_ald <- mock$ald |> filter(!(variable == "AGE" & statistic == "sd"))
  
  expect_error(
    maic.boot(
      ipd = mock$ipd,
      indices = 1:nrow(mock$ipd),
      outcome_model = y ~ trt,
      balance_model = ~ AGE + SCORE,
      family = binomial(),
      ald = bad_ald,
      trt_var = "trt",
      moments = 2,
      int = FALSE
    ),
    "Both 'mean' and 'sd' must be in the ALD to balance the variance of: AGE"
  )
})

# Tests for Interactions (int = TRUE)
test_that("maic.boot successfully runs with int = TRUE when interaction target is provided", {
  mock <- generate_mock_data()
  
  result <- maic.boot(
    ipd = mock$ipd,
    indices = 1:nrow(mock$ipd),
    outcome_model = y ~ trt,
    balance_model = ~ AGE + SCORE,
    family = binomial(),
    ald = mock$ald, # This ALD contains the "AGE:SCORE" mean
    trt_var = "trt",
    moments = 1,
    int = TRUE
  )
  
  expect_true(!is.na(result["ESS"]))
})

test_that("maic.boot fails with int = TRUE if interaction target is missing", {
  mock <- generate_mock_data()
  
  # Remove the interaction term from ALD
  bad_ald <- mock$ald |> filter(variable != "AGE:SCORE")
  
  expect_error(
    maic.boot(
      ipd = mock$ipd,
      indices = 1:nrow(mock$ipd),
      outcome_model = y ~ trt,
      balance_model = ~ AGE + SCORE,
      family = binomial(),
      ald = bad_ald,
      trt_var = "trt",
      moments = 1,
      int = TRUE
    ),
    "Target statistic not found in ALD for covariate: AGE:SCORE"
  )
})

# Test calc_maic passes parameters correctly
test_that("calc_maic processes moments and int arguments via the strategy object", {
  mock <- generate_mock_data()
  
  # Mock strategy object
  mock_strategy <- list(
    n_boot = 5, # keep low for fast testing
    balance_model = ~ AGE + SCORE,
    outcome_model = y ~ trt,
    family = binomial(),
    trt_var = "trt",
    moments = 2,  # Testing pass-through
    int = TRUE    # Testing pass-through
  )
  
  # Mock analysis_params object
  mock_analysis_params <- list(
    ipd = mock$ipd,
    ald = mock$ald,
    verbose = FALSE
  )
  
  # Run the top-level calculation
  results <- calc_maic(
    strategy = mock_strategy, 
    analysis_params = mock_analysis_params
  )
  
  # Assertions
  expect_type(results, "list")
  expect_true(all(c("means", "model") %in% names(results)))
  expect_length(results$model$weights, 100) # Should return weights for all 100 IPD patients
  expect_true(is.numeric(results$model$ESS))
})
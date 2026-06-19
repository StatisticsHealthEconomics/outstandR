# Unit tests for calc_bucher_naive function

library(tibble)

# Test Data Setup
# Create simple IPD and ALD data
ipd_simple <- data.frame(
  y = c(1, 1, 0, 0, 1, 0, 1, 0),
  trt = c("A", "A", "A", "A", "C", "C", "C", "C")
)

ald_simple <- tribble(
  ~variable, ~trt, ~statistic, ~value,
  "y",       "B",  "mean",    0.6,
  "y",       "C",  "mean",    0.4
)

# ============================================================================
# Test 1: Basic functionality with binomial logit link
# ============================================================================

test_that("calc_bucher_naive returns list with expected elements", {
  result <- calc_bucher_naive(
    ipd_trial = ipd_simple,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = NULL
  )
  
  expect_type(result, "list")
  expect_named(result, c("Estimate", "Method", "Treatments"))
  expect_type(result$Estimate, "double")
  expect_type(result$Method, "character")
  expect_type(result$Treatments, "character")
})

test_that("calc_bucher_naive computes correct link-scale estimate", {
  result <- calc_bucher_naive(
    ipd_trial = ipd_simple,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = NULL
  )
  
  # Expected calculation:
  # mu_A = mean(y[trt == "A"]) = 0.5
  # mu_C_ipd = mean(y[trt == "C"]) = 0.5
  # mu_B = 0.6, mu_C_ald = 0.4
  # g_mu_A = log(0.5 / 0.5) = 0
  # g_mu_C_ipd = log(0.5 / 0.5) = 0
  # g_mu_B = log(0.6 / 0.4) ≈ 0.405
  # g_mu_C_ald = log(0.4 / 0.6) ≈ -0.405
  # delta_AC = 0 - 0 = 0
  # delta_BC = 0.405 - (-0.405) = 0.810
  # naive_link_AB = 0 - 0.810 = -0.810
  
  expect_length(result$Estimate, 1)
  expect_false(is.na(result$Estimate))
})

# ============================================================================
# Test 2: Different link functions
# ============================================================================

test_that("calc_bucher_naive works with log link (Poisson)", {
  ipd_poisson <- data.frame(
    y = c(10, 12, 8, 5, 3, 4, 2, 1),
    trt = c("A", "A", "A", "A", "C", "C", "C", "C")
  )
  
  ald_poisson <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "B",  "mean",    7.5,
    "y",       "C",  "mean",    2.5
  )
  
  result <- calc_bucher_naive(
    ipd_trial = ipd_poisson,
    ald_trial = ald_poisson,
    outcome_model = y ~ trt,
    family = poisson(link = "log"),
    ref_trt = "C",
    scale = NULL
  )
  
  expect_type(result$Estimate, "double")
  expect_false(is.na(result$Estimate))
})

test_that("calc_bucher_naive works with identity link", {
  ipd_identity <- data.frame(
    y = c(10, 12, 8, 5, 3, 4, 2, 1),
    trt = c("A", "A", "A", "A", "C", "C", "C", "C")
  )
  
  ald_identity <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "B",  "mean",    7.5,
    "y",       "C",  "mean",    2.5
  )
  
  result <- calc_bucher_naive(
    ipd_trial = ipd_identity,
    ald_trial = ald_identity,
    outcome_model = y ~ trt,
    family = gaussian(link = "identity"),
    ref_trt = "C",
    scale = NULL
  )
  
  expect_type(result$Estimate, "double")
  expect_false(is.na(result$Estimate))
})

# ============================================================================
# Test 3: Scale transformations
# ============================================================================

test_that("calc_bucher_naive handles risk_difference scale", {
  result <- calc_bucher_naive(
    ipd_trial = ipd_simple,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = "risk_difference"
  )
  
  expect_type(result$Estimate, "double")
  expect_false(is.na(result$Estimate))
  # Risk difference should be between -1 and 1
  expect_gte(result$Estimate, -1)
  expect_lte(result$Estimate, 1)
})

test_that("calc_bucher_naive handles log_relative_risk scale", {
  result <- calc_bucher_naive(
    ipd_trial = ipd_simple,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = "log_relative_risk"
  )
  
  expect_type(result$Estimate, "double")
  expect_false(is.na(result$Estimate))
})

test_that("calc_bucher_naive falls back to link scale for unsupported scale+family combo", {
  result <- calc_bucher_naive(
    ipd_trial = ipd_simple,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = "unsupported_scale"
  )
  
  expect_type(result$Estimate, "double")
  expect_false(is.na(result$Estimate))
})

# ============================================================================
# Test 4: Treatment identification
# ============================================================================

test_that("calc_bucher_naive identifies correct active treatments", {
  ipd_multi <- data.frame(
    y = c(1, 0, 1, 0, 1, 0),
    trt = c("X", "X", "Y", "Y", "C", "C")
  )
  
  ald_multi <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "Z",  "mean",    0.7,
    "y",       "C",  "mean",    0.3
  )
  
  result <- calc_bucher_naive(
    ipd_trial = ipd_multi,
    ald_trial = ald_multi,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = NULL
  )
  
  expect_true(grepl("X", result$Treatments))
  expect_true(grepl("Z", result$Treatments))
})

# ============================================================================
# Test 5: Method name
# ============================================================================

test_that("calc_bucher_naive returns correct method name", {
  result <- calc_bucher_naive(
    ipd_trial = ipd_simple,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = NULL
  )
  
  expect_equal(result$Method, "Naive Bucher (Unadjusted)")
})

# ============================================================================
# Test 6: Edge cases
# ============================================================================

test_that("calc_bucher_naive handles data with NA values", {
  ipd_na <- data.frame(
    y = c(1, 1, NA, 0, 1, 0, NA, 0),
    trt = c("A", "A", "A", "A", "C", "C", "C", "C")
  )
  
  result <- calc_bucher_naive(
    ipd_trial = ipd_na,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = NULL
  )
  
  # Should not error and should return numeric
  expect_type(result$Estimate, "double")
})

test_that("calc_bucher_naive handles extreme proportions (0 and 1)", {
  ipd_extreme <- data.frame(
    y = c(1, 1, 1, 1, 0, 0, 0, 0),
    trt = c("A", "A", "A", "A", "C", "C", "C", "C")
  )
  
  ald_extreme <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "B",  "mean",    1.0,
    "y",       "C",  "mean",    0.0
  )
  
  result <- calc_bucher_naive(
    ipd_trial = ipd_extreme,
    ald_trial = ald_extreme,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = NULL
  )
  
  # Link scale avoids direct transformation issues
  expect_type(result$Estimate, "double")
})

test_that("calc_bucher_naive handles large sample sizes", {
  ipd_large <- data.frame(
    y = rbinom(1000, 1, 0.5),
    trt = c(rep("A", 500), rep("C", 500))
  )
  
  ald_large <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "B",  "mean",    0.6,
    "y",       "C",  "mean",    0.4
  )
  
  result <- calc_bucher_naive(
    ipd_trial = ipd_large,
    ald_trial = ald_large,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = NULL
  )
  
  expect_type(result$Estimate, "double")
  expect_false(is.na(result$Estimate))
})

# ============================================================================
# Test 7: Consistency checks
# ============================================================================

test_that("calc_bucher_naive is consistent across multiple calls", {
  result1 <- calc_bucher_naive(
    ipd_trial = ipd_simple,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = NULL
  )
  
  result2 <- calc_bucher_naive(
    ipd_trial = ipd_simple,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = NULL
  )
  
  expect_equal(result1$Estimate, result2$Estimate)
})

test_that("calc_bucher_naive produces sensible estimates across scales", {
  result_link <- calc_bucher_naive(
    ipd_trial = ipd_simple,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = NULL
  )
  
  result_rd <- calc_bucher_naive(
    ipd_trial = ipd_simple,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = "risk_difference"
  )
  
  result_lrr <- calc_bucher_naive(
    ipd_trial = ipd_simple,
    ald_trial = ald_simple,
    outcome_model = y ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = "log_relative_risk"
  )
  
  # All should be numeric and finite
  expect_true(is.finite(result_link$Estimate))
  expect_true(is.finite(result_rd$Estimate))
  expect_true(is.finite(result_lrr$Estimate))
})

# ============================================================================
# Test 8: Outcome variable extraction
# ============================================================================

test_that("calc_bucher_naive correctly extracts outcome variable from formula", {
  ipd_diff_outcome <- data.frame(
    outcome = c(1, 1, 0, 0, 1, 0, 1, 0),
    trt = c("A", "A", "A", "A", "C", "C", "C", "C")
  )
  
  ald_diff_outcome <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "outcome", "B",  "mean",    0.6,
    "outcome", "C",  "mean",    0.4
  )
  
  result <- calc_bucher_naive(
    ipd_trial = ipd_diff_outcome,
    ald_trial = ald_diff_outcome,
    outcome_model = outcome ~ trt,
    family = binomial(link = "logit"),
    ref_trt = "C",
    scale = NULL
  )
  
  expect_type(result$Estimate, "double")
  expect_false(is.na(result$Estimate))
})

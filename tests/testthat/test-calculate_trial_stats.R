# test-calculate_trial_stats

library(tibble)

# variances

test_that("calculate_trial_variance_binary works correctly", {
  ## mock data
  ald <- data.frame(
    variable = c("y", NA, "y", NA),
    trt = c("B", "B", "C", "C"),
    statistic = c("sum", "N", "sum", "N"),
    value = c(30, 100, 20, 100)
  )
  
  # log_odds
  expected <- 1 / 30 + 1 / (100 - 30)
  expect_equal(calculate_trial_variance_binary(ald, "B", "log_odds"),
               expected)
  
  # log_relative_risk
  expected <- 1 / 30 - 1 / 100
  expect_equal(calculate_trial_variance_binary(ald, "B", "log_relative_risk"),
               expected)
  
  # risk_difference
  expected <- 30 * (1 - 30 / 100) / 100
  expect_equal(calculate_trial_variance_binary(ald, "B", "risk_difference"),
               expected)
  
  # delta_z
  expected <- 1 / 30 + 1 / (100 - 30)
  expect_equal(calculate_trial_variance_binary(ald, "B", "delta_z"),
               expected)
  
  # log_relative_risk_rare_events
  expected <- 1 / 30 - 1 / 100
  expect_equal(
    calculate_trial_variance_binary(ald, "B", "log_relative_risk_rare_events"),
    expected
  )
})

test_that("calculate_trial_variance_continuous works correctly", {
  
  ald <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "B",  "mean",     2.5,
    "y",       "B",  "sd",       1.2,
    NA,        "B",  "N",        50
  )
  
  # log_odds
  expected <- pi^2 / 3 * (1 / 50)
  expect_equal(calculate_trial_variance_continuous(ald, "B", "log_odds"),
               expected)
  
  # log_relative_risk
  expected <- log(2.5)
  expect_equal(calculate_trial_variance_continuous(ald, "B", "log_relative_risk"),
               expected)
  
  # risk_difference
  expected <- (1.2^2) / 50
  expect_equal(calculate_trial_variance_continuous(ald, "B", "mean_difference"),
               expected)
})

# means

test_that("calculate_trial_mean_binary works correctly", {

  ald <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "B",  "sum",     30,
    NA,        "B",  "N",       100
  )
  
  p <- 30 / 100
  
  # log_odds
  expected <- qlogis(p)
  expect_equal(calculate_trial_mean_binary(ald, "B", "log_odds"), expected)
  
  # risk_difference
  expected <- p
  expect_equal(calculate_trial_mean_binary(ald, "B", "risk_difference"),
               expected)
  
  # delta_z
  expected <- qnorm(p)
  expect_equal(calculate_trial_mean_binary(ald, "B", "delta_z"), expected)
  
  # log_relative_risk_rare_events
  expected <- log(-log(1 - p))
  expect_equal(
    calculate_trial_mean_binary(ald, "B", "log_relative_risk_rare_events"),
    expected
  )
  
  # log_relative_risk
  expected <- log(p)
  expect_equal(calculate_trial_mean_binary(ald, "B", "log_relative_risk"),
               expected)
})



test_that("calculate_trial_mean_continuous works correctly", {
  
  ald <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "B",  "mean",     2.5,
    "y",       "B",  "sd",       1.2,
    NA,        "B",  "N",        50
  )

  # log_odds
  expected <- log(2.5)
  expect_equal(calculate_trial_mean_continuous(ald, "B", "log_odds"),
               expected)
  
  # risk_difference
  expected <- 2.5
  expect_equal(calculate_trial_mean_continuous(ald, "B", "mean_difference"),
               expected)
  
  # delta_z
  expected <- 2.5 / 1.2
  expect_equal(calculate_trial_mean_continuous(ald, "B", "delta_z"),
               expected)
})



test_that("marginal_variance works correctly", {

  ald <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "B",  "sum",     30,
    "y",       "C",  "sum",     40,
    NA,        "B",  "N",       100,
    NA,        "C",  "N",       120
  )

  expected <- sum(
    calculate_trial_variance_binary(ald, "B", "log_odds"),
    calculate_trial_variance_binary(ald, "C", "log_odds")
  )
  
  expect_equal(
    marginal_variance(ald, 
                      ref_trt = "C", comp_trt = "B", 
                      scale = "log_odds", family = "binomial"),
    expected)
})


test_that("marginal_treatment_effect works correctly", {
  
  ald <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "B",  "sum",     30,
    "y",       "C",  "sum",     40,
    NA,        "B",  "N",       100,
    NA,        "C",  "N",       120
  )
  
  meanB <- calculate_trial_mean_binary(ald, "B", "log_odds")
  meanC <- calculate_trial_mean_binary(ald, "C", "log_odds")
  
  expected <- meanB - meanC 
  
  expect_equal(
    marginal_treatment_effect(ald, ref_trt = "C", comp_trt = "B", 
                              scale = "log_odds", family = "binomial"),
    expected)
})

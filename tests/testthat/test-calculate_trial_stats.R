# test-calculate_trial_stats

# variances

test_that("calculate_trial_variance_binary works correctly", {
  ald <- list(
    "y.B.sum" = 30,
    "N.B" = 100,
    "y.C.sum" = 40,
    "N.C" = 120
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
  ald <- list("y.B.bar" = 2.5,
              "y.B.sd" = 1.2,
              "N.B" = 50)
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
  expect_equal(calculate_trial_variance_continuous(ald, "B", "risk_difference"),
               expected)
})

# means

test_that("calculate_trial_mean_binary works correctly", {
  ald <- list("y.B.sum" = 30,
              "N.B" = 100)
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
  ald <- list("y.B.bar" = 2.5,
              "y.B.sd" = 1.2,
              "N.B" = 50)
  
  # log_odds
  expected <- log(2.5)
  expect_equal(calculate_trial_mean_continuous(ald, "B", "log_odds"),
               expected)
  
  # risk_difference
  expected <- 2.5
  expect_equal(calculate_trial_mean_continuous(ald, "B", "risk_difference"),
               expected)
  
  # delta_z
  expected <- 2.5 / 1.2
  expect_equal(calculate_trial_mean_continuous(ald, "B", "delta_z"),
               expected)
})



test_that("marginal_variance works correctly", {
  ald <- list(
    "y.B.sum" = 30,
    "N.B" = 100,
    "y.C.sum" = 40,
    "N.C" = 120
  )
  treatments <- list("B", "C")
  
  expected <- sum(
    calculate_trial_variance_binary(ald, "B", "log_odds"),
    calculate_trial_variance_binary(ald, "C", "log_odds")
  )
  
  expect_equal(marginal_variance(ald, treatments, "log_odds", "binomial"),
               expected)
})


test_that("marginal_treatment_effect works correctly", {
  ald <- list(
    "y.B.sum" = 30,
    "N.B" = 100,
    "y.C.sum" = 40,
    "N.C" = 120
  )
  treatments <- list("B", "C")
  
  meanB <- calculate_trial_mean_binary(ald, "B", "log_odds")
  meanC <- calculate_trial_mean_binary(ald, "C", "log_odds")
  
  expected <- meanB - meanC 
  expect_equal(marginal_treatment_effect(ald, treatments, "log_odds", "binomial"),
               expected)
})

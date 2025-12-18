# This test suite checks three key scenarios:
#   
# Auto- calculation: Providing just the distribution name (e.g., "gamma") and
#                    verifying that outstandR correctly converts the ALD's Mean/SD into the correct Shape/Rate.
#
# Manual override: Providing both distribution and parameters to ensure the user can still force specific values.
#
# Integration: Ensuring these arguments work when passed through the top-level outstandR() function.

# test simulate_ALD_pseudo_pop()

test_that("Method of Moments Auto-Calculation (Gamma & Beta)", {
  
  # --- Setup Mock Data ---
  # Create a mock ALD with known Mean/SD
  # Variable 1: Cost (Gamma-like: Mean=1000, SD=500)
  # Variable 2: Utility (Beta-like: Mean=0.7, SD=0.1)
  ald_mock <- data.frame(
    variable = c("cost", "cost", "util", "util"),
    statistic = c("mean", "sd", "mean", "sd"),
    value = c(1000, 500, 0.7, 0.1),
    trt = NA
  )
  
  # Formula (Treatment is not in ALD covariate rows)
  form <- outcome ~ cost + util*trt + trt
  trt_var_name <- "trt"
  
  # --- TEST 1: Gamma Auto-Conversion ---
  
  # Expected Method of Moments for Gamma:
  # shape = mean^2 / sd^2 = 1000^2 / 500^2 = 4
  # rate = mean / sd^2 = 1000 / 250000 = 0.004
  # Expected Mean = shape/rate = 1000
  # Expected SD = sqrt(shape/rate^2) = 500
  
  set.seed(123)
  
  rho_mat <- diag(2)  # independence
  rho_names <- c("cost", "util")
  dimnames(rho_mat) <- list(rho_names, rho_names)
  
  sim_gamma <- simulate_ALD_pseudo_pop(
    formula = form,
    ald = ald_mock,
    trt_var = trt_var_name,
    rho = rho_mat,
    N = 10000,
    marginal_distns = c(cost = "gamma", util = "norm") # 'util' default norm
  )
  
  # simulated data matches the ALD moments (allow small simulation error)
  expect_equal(mean(sim_gamma$cost), 1000, tolerance = 0.05)
  expect_equal(sd(sim_gamma$cost), 500, tolerance = 0.05)
  
  # NOT Gaussian (skewness check or range check)
  # Gamma(4, 0.004) should be strictly positive
  expect_true(min(sim_gamma$cost) > 0)
  
  
  # --- TEST 2: Beta Auto-Conversion ---
  
  # Expected Method of Moments for Beta:
  # term = (mean * (1-mean) / var) - 1
  #      = (0.7 * 0.3 / 0.01) - 1 = (0.21 / 0.01) - 1 = 21 - 1 = 20
  # shape1 = mean * term = 0.7 * 20 = 14
  # shape2 = (1-mean) * term = 0.3 * 20 = 6
  
  set.seed(123)
  
  sim_beta <- simulate_ALD_pseudo_pop(
    formula = form,
    ald = ald_mock,
    trt_var = trt_var_name,
    rho = rho_mat,
    N = 10000,
    marginal_distns = c(cost = "norm", util = "beta")
  )
  
  # Check moments
  expect_equal(mean(sim_beta$util), 0.7, tolerance = 0.01)
  expect_equal(sd(sim_beta$util), 0.1, tolerance = 0.01)
  
  # Verify bounded [0,1]
  expect_true(max(sim_beta$util) <= 1)
  expect_true(min(sim_beta$util) >= 0)
})


test_that("Manual Override of ALD Parameters", {
  
  ald_mock <- data.frame(
    variable = c("age", "age"),
    statistic = c("mean", "sd"),
    value = c(50, 10),  # ALD Mean=50
    trt = NA
  )
  
  form <- outcome ~ age + trt
  
  # We want to force Mean=70 (Sensitivity Analysis)
  # target_mean=70, target_sd=10 -> shape=49, rate=0.7
  forced_params <- list(
    age = list(shape = 49, rate = 0.7)
  )
  
  set.seed(123)
  
  sim_manual <- simulate_ALD_pseudo_pop(
    formula = form,
    ald = ald_mock,
    trt_var = "trt",
    rho = 1,
    N = 10000,
    marginal_distns = c(age = "gamma"),
    marginal_params = forced_params # <--- OVERRIDE
  )
  
  # Should match Manual (70), NOT ALD (50)
  expect_equal(mean(sim_manual$age), 70, tolerance = 0.5)
  expect_false(abs(mean(sim_manual$age) - 50) < 1)
})


test_that("Integration: Arguments pass through outstandR() wrapper", {
  
  # Load package data
  # Use AC_IPD_binY_contX / BC_ALD_binY_contX (loaded in helper or available in pkg)
  # Variables: EM_cont_1 (mean~0.65, sd~0.39)
  
  # 1. Define strategy with custom distribution
  my_strategy <- strategy_gcomp_ml(
    formula = y ~ PF_cont_1 + PF_cont_2 + trt + trt:EM_cont_1 + trt:EM_cont_2,
    family = binomial(link = "logit"),
    marginal_distns = c(EM_cont_1 = "gamma", 
                        EM_cont_2 = "norm", 
                        PF_cont_1 = "norm", 
                        PF_cont_2 = "norm"),
    N = 200 # Small N for speed
  )
  
  # 2. Run wrapper
  # If the arguments weren't threaded correctly in R/gcomp_ml.R, 
  # this would crash (Gamma needs shape/rate, defaults to norm if missing) 
  # or silently fail to produce Gamma data.
  
  res <- outstandR(
    ipd_trial = AC_IPD_binY_contX,
    ald_trial = BC_ALD_binY_contX,
    strategy = my_strategy
  )
  
  expect_s3_class(res, "outstandR")
  expect_true(is.numeric(res$contrasts$mean))
})

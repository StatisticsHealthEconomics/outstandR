test_that("Sandwich variance estimator works for STC", {
  # Setup data
  # Use binary outcome data available in the package
  ipd <- AC_IPD_binY_contX
  ald <- BC_ALD_binY_contX
  
  # Define Strategy
  form <- y ~ PF_cont_1 + PF_cont_2 + trt + trt:EM_cont_1 + trt:EM_cont_2
  strat <- strategy_stc(formula = form, 
                        family = binomial(link = "logit"), 
                        trt_var = "trt")
  
  params <- list(ipd = ipd, ald = ald, 
                 ref_trt = "C", ipd_comp = "A", 
                 scale = "log_odds")
  
  # 1. Run with Naive Variance (Default)
  res_naive <- calc_IPD_stats(strat, params, var_method = "sample")
  
  # 2. Run with Sandwich Variance
  res_sandwich <- calc_IPD_stats(strat, params, var_method = "sandwich")
  
  # Checks
  expect_true(is.numeric(res_sandwich$contrasts$var))
  expect_gt(res_sandwich$contrasts$var, 0)
  
  # Sandwich variance should essentially never equal naive sample variance exactly
  expect_false(res_sandwich$contrasts$var == res_naive$contrasts$var)
  
  # Check structure matches
  expect_named(res_sandwich$contrasts, c("mean", "var"))
})

test_that("Sandwich variance estimator works for MAIC", {
  ipd <- AC_IPD_binY_contX
  ald <- BC_ALD_binY_contX
  
  # Define Strategy (MAIC)
  form <- y ~ PF_cont_1 + PF_cont_2 + trt + trt:EM_cont_1 + trt:EM_cont_2
  strat <- strategy_maic(formula = form, 
                         family = binomial(link = "logit"), 
                         trt_var = "trt")
  
  params <- list(ipd = ipd, ald = ald, 
                 ref_trt = "C", ipd_comp = "A", 
                 scale = "risk_difference") # Test different scale
  
  # Run Sandwich
  res_sandwich <- calc_IPD_stats(strat, params, var_method = "sandwich")
  
  expect_true(is.numeric(res_sandwich$contrasts$var))
  expect_gt(res_sandwich$contrasts$var, 0)
  
  # Verify that scale transformation (Delta Method) worked
  # Risk difference should have different magnitude variance than log odds
  expect_true(res_sandwich$contrasts$var < 1) # RD variance is usually small < 1
})

test_that("Sandwich variance estimator works for G-Computation (ML)", {
  ipd <- AC_IPD_binY_contX
  ald <- BC_ALD_binY_contX
  
  form <- y ~ PF_cont_1 + PF_cont_2 + trt + trt:EM_cont_1 + trt:EM_cont_2
  strat <- strategy_gcomp_ml(formula = form, 
                             family = binomial(link = "logit"), 
                             trt_var = "trt",
                             N = 100) # Small N for speed
  
  params <- list(ipd = ipd, ald = ald, 
                 ref_trt = "C", ipd_comp = "A", 
                 scale = "log_relative_risk")
  
  res_sandwich <- calc_IPD_stats(strat, params, var_method = "sandwich")
  
  expect_true(is.numeric(res_sandwich$contrasts$var))
  expect_gt(res_sandwich$contrasts$var, 0)
})

test_that("Sandwich estimator throws error for unsupported methods", {
  # Bayesian G-comp does not support sandwich (uses posterior variance)
  strat_bayes <- strategy_gcomp_bayes(formula = y ~ trt, 
                                      family = gaussian(), 
                                      trt_var = "trt", N = 10)
  
  params <- list(ipd = AC_IPD_contY_mixedX, ald = BC_ALD_contY_mixedX,
                 ref_trt = "C", ipd_comp = "A", scale = "mean_difference")
  
  expect_error(
    calc_IPD_stats(strat_bayes, params, var_method = "sandwich"),
    regexp = "Sandwich variance not implemented"
  )
})

test_that("Internal helper: Numerical Gradient works", {
  # Simple quadratic function y = x^2
  # Gradient at x=2 should be 2x = 4
  func <- function(x) x^2
  grad <- outstandR:::num_grad(func, x = 2)
  expect_equal(grad, 4, tolerance = 1e-4)
  
  # Multivariate: y = x1^2 + 3*x2
  # Grad = [2*x1, 3]
  func_multi <- function(x) x[1]^2 + 3*x[2]
  grad_multi <- outstandR:::num_grad(func_multi, x = c(2, 5))
  expect_equal(grad_multi, c(4, 3), tolerance = 1e-4)
})
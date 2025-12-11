test_that("Sandwich variance estimator works for STC", {
  ipd <- AC_IPD_binY_contX
  ald <- BC_ALD_binY_contX
  
  form <- y ~ PF_cont_1 + PF_cont_2 + trt + trt:EM_cont_1 + trt:EM_cont_2
  strat <- strategy_stc(formula = form, 
                        family = binomial(link = "logit"), 
                        trt_var = "trt")
  
  # Base parameters
  params <- list(ipd = ipd, ald = ald, 
                 ref_trt = "C", ipd_comp = "A", 
                 scale = "log_odds")
  
  # 1. Run with Naive Variance (Default)
  res_naive <- calc_IPD_stats(strat, params) # Defaults to "sample" internally
  
  # 2. Run with Sandwich Variance
  params_sand <- params
  params_sand$var_method <- "sandwich"
  res_sandwich <- calc_IPD_stats(strat, params_sand)
  
  # Checks
  expect_true(is.numeric(res_sandwich$contrasts$var))
  expect_gt(res_sandwich$contrasts$var, 0)
  expect_false(res_sandwich$contrasts$var == res_naive$contrasts$var)
})

test_that("Sandwich variance estimator works for MAIC", {
  ipd <- AC_IPD_binY_contX
  ald <- BC_ALD_binY_contX
  
  form <- y ~ PF_cont_1 + PF_cont_2 + trt + trt:EM_cont_1 + trt:EM_cont_2
  strat <- strategy_maic(formula = form, 
                         family = binomial(link = "logit"), 
                         trt_var = "trt")
  
  # Specify method inside params
  params <- list(ipd = ipd, ald = ald, 
                 ref_trt = "C", ipd_comp = "A", 
                 scale = "risk_difference",
                 var_method = "sandwich")
  
  res_sandwich <- calc_IPD_stats(strat, params)
  
  expect_true(is.numeric(res_sandwich$contrasts$var))
  expect_gt(res_sandwich$contrasts$var, 0)
})

test_that("Sandwich variance works via top-level outstandR() wrapper", {
  # This tests that the wrapper correctly passes the argument down into params
  
  # 1. Naive run
  res_naive <- outstandR(
    ipd_trial = AC_IPD_binY_contX,
    ald_trial = BC_ALD_binY_contX,
    strategy = strategy_stc(formula = y ~ PF_cont_1 + trt, family = binomial()),
    var_method = "sample"
  )
  
  # 2. Robust run
  res_robust <- outstandR(
    ipd_trial = AC_IPD_binY_contX,
    ald_trial = BC_ALD_binY_contX,
    strategy = strategy_stc(formula = y ~ PF_cont_1 + trt, family = binomial()),
    var_method = "sandwich" 
  )
  
  expect_false(res_naive$contrasts$var == res_robust$contrasts$var)
  expect_equal(res_robust$scale, "log_odds") # Check other attributes preserved
})
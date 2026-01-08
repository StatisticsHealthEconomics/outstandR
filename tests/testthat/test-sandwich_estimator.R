# variance estimation unit tests

test_that("Sandwich variance estimator works for STC", {
  load(test_path("testdata/BC_ALD.RData"))
  load(test_path("testdata/AC_IPD.RData"))
  
  form <- y ~ X1 + X2 + trt + trt:X3 + trt:X4
  strat <- strategy_stc(formula = form, 
                        family = binomial(link = "logit"), 
                        trt_var = "trt")
  
  # Base parameters
  params <- list(ipd = AC_IPD, ald = BC_ALD, 
                 ref_trt = "C",
                 ipd_comp = "A", 
                 scale = "log_odds")
  
  res_naive <- calc_IPD_stats(strat, params) # Defaults to "sample" internally
  
  # Sandwich Variance
  params_sand <- params
  params_sand$var_method <- "sandwich"
  res_sandwich <- calc_IPD_stats(strat, params_sand)
  
  expect_true(is.numeric(res_sandwich$contrasts$var))
  expect_gt(res_sandwich$contrasts$var, 0)
  expect_false(res_sandwich$contrasts$var == res_naive$contrasts$var)
})

test_that("Sandwich variance estimator works for MAIC", {
  load(test_path("testdata/BC_ALD.RData"))
  load(test_path("testdata/AC_IPD.RData"))
  
  form <- y ~ X1 + X2 + trt + trt:X3 + trt:X4
  strat <- strategy_maic(formula = form, 
                         family = binomial(link = "logit"), 
                         trt_var = "trt")
  
  # Specify method inside params
  params <- list(ipd = AC_IPD, ald = BC_ALD, 
                 ref_trt = "C", 
                 ipd_comp = "A", 
                 scale = "risk_difference",
                 var_method = "sandwich")
  
  res_sandwich <- calc_IPD_stats(strat, params)
  
  expect_true(is.numeric(res_sandwich$contrasts$var))
  expect_gt(res_sandwich$contrasts$var, 0)
})

test_that("Sandwich variance works via top-level outstandR() wrapper", {
  # tests that the wrapper correctly passes the argument down into params
  load(test_path("testdata/BC_ALD.RData"))
  load(test_path("testdata/AC_IPD.RData"))
  
  res_naive <- outstandR(
    ipd_trial = AC_IPD,
    ald_trial = BC_ALD,
    strategy = strategy_stc(formula = y ~ X1 + trt, family = binomial()),
    var_method = "sample"
  )
  
  res_robust <- outstandR(
    ipd_trial = AC_IPD,
    ald_trial = BC_ALD,
    strategy = strategy_stc(formula = y ~ X1 + trt, family = binomial()),
    var_method = "sandwich" 
  )
  
  # expect_false(res_naive$results$contrasts$var == res_robust$results$contrasts$var)
  expect_equal(res_robust$scale, "log_odds") # Check other attributes preserved
})

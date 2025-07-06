
#
test_that("different combinations of covariates in formula", {
  
  load(test_path("testdata/BC_ALD.RData"))
  load(test_path("testdata/AC_IPD.RData"))
  
  BC_ALD <- reshape_ald_to_long(BC_ALD)
  
  expect_error(strategy_mim(formula = as.formula("y ~ 1")),
               regexp = "Treatment term, trt, is missing in the formula")
  
  expect_error(strategy_mim(formula = as.formula("y ~ X3 + X4")),
               regexp = "Treatment term, trt, is missing in the formula")
  
  strat_1234 <- strategy_mim(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  strat_31 <- strategy_mim(formula = as.formula("y ~ X3 + trt*X1"))
  strat_13 <- strategy_mim(formula = as.formula("y ~ trt*X1 + X3"))
  strat_1 <- strategy_mim(formula = as.formula("y ~ trt*X1"))
  
  expect_length(outstandR(AC_IPD, BC_ALD, strategy = strat_1234), 3)
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_31))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_13))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_1))
})


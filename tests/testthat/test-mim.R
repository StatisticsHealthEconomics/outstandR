
#
test_that("different combinations of covariates in formula", {
  
  load(test_path("testdata/BC_ALD.RData"))
  load(test_path("testdata/AC_IPD.RData"))
  
  BC_ALD <- reshape_ald_to_long(BC_ALD)
  AC_IPD$trt <- factor(AC_IPD$trt, labels = c("C", "A"))  # from 0, 1
  
  expect_error(strategy_mim(formula = as.formula("y ~ 1")),
               regexp = "Treatment term 'trt' is missing in the formula")
  
  expect_message(strategy_mim(formula = as.formula("y ~ X3 + X4")),
               regexp = "Treatment is guessed as:")
  
  strat_1234 <- strategy_mim(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  strat_31 <- strategy_mim(formula = as.formula("y ~ X3 + trt*X1"))
  strat_13 <- strategy_mim(formula = as.formula("y ~ trt*X1 + X3"))
  strat_1 <- strategy_mim(formula = as.formula("y ~ trt*X1"))
  
  expect_length(outstandR(AC_IPD, BC_ALD, strategy = strat_1234), 2)
  expect_named(outstandR(AC_IPD, BC_ALD, strategy = strat_1234),
               expected = c("contrasts", "absolute"))
})


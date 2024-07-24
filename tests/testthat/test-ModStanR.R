#


test_that("ModStanR basics", {
  expect_equal(2 * 2, 4)
})


test_that("ModStanR errors", {
  expect_equal(2 * 2, 4)
})


#
test_that("different combinations of covariates in formula", {
  
  load(test_path("testdata/BC_ALD.RData"))
  load(test_path("testdata/AC_IPD.RData"))
  
  # maic
  expect_error(strategy_maic(formula = as.formula("y ~ 1")),
               regexp = "Treatment term, trt, is missing in the formula")

  expect_error(strategy_maic(formula = as.formula("y ~ X3 + X4")),
               regexp = "Treatment term, trt, is missing in the formula")
  
  strat_1234 <- strategy_maic(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  strat_31 <- strategy_maic(formula = as.formula("y ~ X3 + trt*X1"))
  strat_13 <- strategy_maic(formula = as.formula("y ~ trt*X1 + X3"))
  strat_1 <- strategy_maic(formula = as.formula("y ~ trt*X1"))
  
  expect_length(ModStanR(AC_IPD, BC_ALD, strategy = strat_1234), 3)
  expect_equal(ModStanR(AC_IPD, BC_ALD, strategy = strat_31))
  expect_equal(ModStanR(AC_IPD, BC_ALD, strategy = strat_13))
  expect_equal(ModStanR(AC_IPD, BC_ALD, strategy = strat_1))

  # stc
  expect_error(strategy_stc(formula = as.formula("y ~ 1")),
               regexp = "Treatment term, trt, is missing in the formula")
  
  expect_error(strategy_stc(formula = as.formula("y ~ X3 + X4")),
               regexp = "Treatment term, trt, is missing in the formula")
  
  strat_1234 <- strategy_stc(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  strat_31 <- strategy_stc(formula = as.formula("y ~ X3 + trt*X1"))
  strat_13 <- strategy_stc(formula = as.formula("y ~ trt*X1 + X3"))
  strat_1 <- strategy_stc(formula = as.formula("y ~ trt*X1"))

  expect_equal(ModStanR(AC_IPD, BC_ALD, strategy = strat_1234))
  expect_equal(ModStanR(AC_IPD, BC_ALD, strategy = strat_31))
  expect_equal(ModStanR(AC_IPD, BC_ALD, strategy = strat_13))
  expect_equal(ModStanR(AC_IPD, BC_ALD, strategy = strat_1))
})

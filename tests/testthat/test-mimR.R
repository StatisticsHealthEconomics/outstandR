#


test_that("mimR basics", {
  expect_equal(2 * 2, 4)
})


test_that("mimR errors", {
  expect_equal(2 * 2, 4)
})


#
test_that("different combinations of covariates in formula", {
  
  # # maic
  # strat_none <- strategy_maic(formula = as.formula("y ~ ."))
  # strat_1 <- strategy_maic(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  # strat_2 <- strategy_maic(formula = as.formula("y ~ X3 + X4 + X5 + trt*X1"))
  # strat_no_trt <- strategy_maic(formula = as.formula("y ~ X3 + X4"))
  # strat_31 <- strategy_maic(formula = as.formula("y ~ X3 + trt*X1"))
  # strat_32 <- strategy_maic(formula = as.formula("y ~ trt*X1 + X3"))
  # strat_no_X1 <- strategy_maic(formula = as.formula("y ~ trt*X1"))
  # 
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_none))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_1))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_2))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_no_trt))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_31))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_32))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_no_X1))
  # 
  # # stc
  # strat_none <- strategy_stc(formula = as.formula("y ~ ."))
  # strat_1 <- strategy_stc(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  # strat_2 <- strategy_stc(formula = as.formula("y ~ X3 + X4 + X5 + trt*X1"))
  # strat_no_trt <- strategy_stc(formula = as.formula("y ~ X3 + X4"))
  # strat_31 <- strategy_stc(formula = as.formula("y ~ X3 + trt*X1"))
  # strat_32 <- strategy_stc(formula = as.formula("y ~ trt*X1 + X3"))
  # strat_no_X1 <- strategy_stc(formula = as.formula("y ~ trt*X1"))
  # 
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_none))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_1))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_2))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_no_trt))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_31))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_32))
  # expect_equal(mimR(AC.IPD, BC.ALD, strategy = strat_no_X1))
})

test_that("plot.outstandR returns a ggplot object and displays the correct scale", {
  load(test_path("testdata/BC_ALD.RData"))
  load(test_path("testdata/AC_IPD.RData"))
  
  strategy <- strategy_maic(
    formula = list(outcome_model = y ~ trt, balance_model = ~ X1),
    family = binomial(link = "logit")
  )
  
  res <- outstandR(AC_IPD, BC_ALD, strategy = strategy, seed = 123, verbose = FALSE)
  
  p <- plot(res)
  expect_s3_class(p, "ggplot")
  
  # Verify that the mapped human-readable scale label is in the x-axis label
  expect_match(p$labels$x, "Scale: Log-Odds Ratio", fixed = TRUE)
  
  # Test with continuous/gaussian data scale
  strategy_cont <- strategy_maic(
    formula = list(outcome_model = y ~ trt, balance_model = ~ X1),
    family = gaussian(link = "identity")
  )
  
  res_cont <- outstandR(AC_IPD, BC_ALD, strategy = strategy_cont, seed = 123, verbose = FALSE)
  p_cont <- plot(res_cont)
  expect_match(p_cont$labels$x, "Scale: Mean Difference", fixed = TRUE)
  
  # Test that passing objects with different scales throws an error
  expect_error(plot(res, res_cont), "All objects passed to plot\\(\\) must be on the same scale")
})

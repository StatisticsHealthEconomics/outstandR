## originally from chatgpt


## mock data
ald <- data.frame(
  variable = c("y", NA, "y", NA),
  trt = c("B", "B", "C", "C"),
  statistic = c("sum", "N", "sum", "N"),
  value = c(30, 100, 20, 100)
)

strategy <- list(family = binomial(link = "logit"))

## calc_ALD_stats() ----
test_that("calc_ALD_stats() returns mean and variance", {

  res <- calc_ALD_stats(
    strategy,
    list(ald = ald, 
         ref_trt = "C", 
         ald_comp = "B", 
         scale = "log_odds"))
  
  expect_type(res, "list")
  expect_named(res, c("mean", "var"))
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

## marginal_variance() ----
test_that("marginal_variance() calculates correct sum of variances", {
  res <- marginal_variance(
    ald, ref_trt = "C", comp_trt = "B",
    scale = "log_odds",
    family = strategy$family$family)
  
  expect_type(res, "double")
  expect_gt(res, 0)
})

## marginal_treatment_effect() ----
test_that("marginal_treatment_effect() calculates correct difference", {
  res <- marginal_treatment_effect(
    ald, ref_trt = "C", comp_trt = "B",
    scale = "log_odds",
    family = strategy$family$family)
  
  expect_type(res, "double")
})

## Edge Cases for calc_ALD_stats() ----
test_that("calc_ALD_stats() handles NULL or empty ald", {
  expect_error(calc_ALD_stats(strategy, NULL))
  expect_error(calc_ALD_stats(strategy, list()))
})

test_that("calc_ALD_stats() handles missing treatment labels", {
  ald_missing <- data.frame(
    variable = c("y", NA, "y", NA),
    trt = c("B", "B", NA, "C"),
    statistic = c("sum", "N", "sum", "N"),
    value = c(30, 100, 20, 100)
    # C is missing
  )
  
  expect_error(calc_ALD_stats(strategy, ald_missing))
})

test_that("calc_ALD_stats() handles incorrect data types", {
  ald_wrong <- data.frame(
      variable = c("y", NA, "y", NA),
      trt = c("B", "B", "C", "C"),
      statistic = c("sum", "N", "sum", "N"),
      value = c("30", "100", "20", "100")
    )
  
  expect_error(calc_ALD_stats(strategy, ald_wrong))
})

test_that("calc_ALD_stats() handles extreme values", {
  ald_extreme <- data.frame(
      variable = c("y", NA, "y", NA),
      trt = c("B", "B", "C", "C"),
      statistic = c("sum", "N", "sum", "N"),
      value = c(0, 100, 100, 100) # zero and all events
    )

  res <- calc_ALD_stats(
    strategy,
    list(ald = ald_extreme, 
         ref_trt = "C", 
         ald_comp = "B", 
         scale = "log_odds"))
  
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})


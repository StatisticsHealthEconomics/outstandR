## from chatgpt


## Mock Data
ald <- list(
  y.B.sum = 30,
  N.B = 100,
  y.C.sum = 20,
  N.C = 100
)

strategy <- list(family = binomial(link = "logit"))

## calc_ALD_stats() ----
test_that("calc_ALD_stats() returns mean and variance", {
  res <- calc_ALD_stats(strategy, ald, scale = "log_odds")
  expect_type(res, "list")
  expect_named(res, c("mean", "var"))
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

## marginal_variance() ----
test_that("marginal_variance() calculates correct sum of variances", {
  res <- marginal_variance(ald, family = strategy$family$family, scale = "log_odds")
  expect_type(res, "double")
  expect_gt(res, 0)
})

## marginal_treatment_effect() ----
test_that("marginal_treatment_effect() calculates correct difference", {
  res <- marginal_treatment_effect(ald, family = strategy$family$family, scale = "log_odds")
  expect_type(res, "double")
})

## Edge Cases for calc_ALD_stats() ----
test_that("calc_ALD_stats() handles NULL or empty ald", {
  expect_error(calc_ALD_stats(strategy, NULL))
  expect_error(calc_ALD_stats(strategy, list()))
})

test_that("calc_ALD_stats() handles missing treatment labels", {
  ald_missing <- list(
    y.B.sum = 30,
    N.B = 100
    # C is missing
  )
  expect_error(calc_ALD_stats(strategy, ald_missing))
})

test_that("calc_ALD_stats() handles incorrect data types", {
  ald_wrong <- list(
    y.B.sum = "thirty",
    N.B = 100,
    y.C.sum = 20,
    N.C = "one hundred"
  )
  expect_error(calc_ALD_stats(strategy, ald_wrong))
})

test_that("calc_ALD_stats() handles extreme values", {
  ald_extreme <- list(
    y.B.sum = 0,   # Zero events
    N.B = 100,
    y.C.sum = 100, # All events
    N.C = 100
  )
  res <- calc_ALD_stats(strategy, ald_extreme, scale = "log_odds")
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

## Edge Cases for marginal_variance() ----
## one group everyone and the other no-one
test_that("marginal_variance() handles zero and full counts", {
  
  ald_extreme <- list(
    y.B.sum = 0,
    N.B = 100,
    y.C.sum = 100,
    N.C = 100
  )
  
  res <- calc_ALD_stats(strategy, ald_extreme, scale = "log_odds") |> unlist()
  
  expect_true(all(is.finite(res)))
})



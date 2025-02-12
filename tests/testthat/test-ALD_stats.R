## from chatgpt


## Mock Data
ald <- list(
  y.B.sum = 30,
  N.B = 100,
  y.C.sum = 20,
  N.C = 100
)

strategy <- list(family = binomial(link = "logit"))

## ALD_stats() ----
test_that("ALD_stats() returns mean and variance", {
  res <- ALD_stats(strategy, ald)
  expect_type(res, "list")
  expect_named(res, c("mean", "var"))
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

## marginal_variance() ----
test_that("marginal_variance() calculates correct sum of variances", {
  res <- marginal_variance(ald, family = strategy$family)
  expect_type(res, "double")
  expect_gt(res, 0)
})

## marginal_treatment_effect() ----
test_that("marginal_treatment_effect() calculates correct difference", {
  res <- marginal_treatment_effect(ald, family = strategy$family)
  expect_type(res, "double")
})

## trial_variance() ----
test_that("trial_variance() returns correct variance for treatment", {
  res_B <- trial_variance(ald, "B", strategy$family)
  res_C <- trial_variance(ald, "C", strategy$family)
  expect_type(res_B, "double")
  expect_type(res_C, "double")
  expect_gt(res_B, 0)
  expect_gt(res_C, 0)
})

## trial_treatment_effect() ----
test_that("trial_treatment_effect() returns correct log-odds ratio", {
  res <- trial_treatment_effect(ald, "B", strategy$family)
  expect_type(res, "double")
})

## link_transform() ----
test_that("link_transform() works with logit link", {
  res <- link_transform(0.3, strategy$family)
  expect_type(res, "double")
})

## link_transform_var() ----
test_that("link_transform_var() returns variance for logit link", {
  res <- link_transform_var(30, 100, strategy$family)
  expect_type(res, "double")
  expect_gt(res, 0)
})

## Edge Cases for ALD_stats() ----
test_that("ALD_stats() handles NULL or empty ald", {
  expect_error(ALD_stats(strategy, NULL))
  expect_error(ALD_stats(strategy, list()))
})

test_that("ALD_stats() handles missing treatment labels", {
  ald_missing <- list(
    y.B.sum = 30,
    N.B = 100
    # C is missing
  )
  expect_error(ALD_stats(strategy, ald_missing))
})

test_that("ALD_stats() handles incorrect data types", {
  ald_wrong <- list(
    y.B.sum = "thirty",
    N.B = 100,
    y.C.sum = 20,
    N.C = "one hundred"
  )
  expect_error(ALD_stats(strategy, ald_wrong))
})

test_that("ALD_stats() handles extreme values", {
  ald_extreme <- list(
    y.B.sum = 0,   # Zero events
    N.B = 100,
    y.C.sum = 100, # All events
    N.C = 100
  )
  res <- ALD_stats(strategy, ald_extreme)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

## Edge Cases for marginal_variance() ----
test_that("marginal_variance() handles zero and full counts", {
  ald_extreme <- list(
    y.B.sum = 0,
    N.B = 100,
    y.C.sum = 100,
    N.C = 100
  )
  res <- marginal_variance(ald_extreme, family = strategy$family)
  expect_true(is.finite(res))
})

## Edge Cases for marginal_treatment_effect() ----
test_that("marginal_treatment_effect() handles edge case differences", {
  ald_extreme <- list(
    y.B.sum = 100,
    N.B = 100,
    y.C.sum = 0,
    N.C = 100
  )
  res <- marginal_treatment_effect(ald_extreme, family = strategy$family)
  expect_true(is.finite(res))
})

## Edge Cases for link_transform() ----
test_that("link_transform() handles different link functions", {
  family_log <- binomial(link = "log")
  family_identity <- binomial(link = "identity")
  res_log <- link_transform(0.3, family_log)
  res_identity <- link_transform(0.3, family_identity)
  
  expect_type(res_log, "double")
  expect_type(res_identity, "double")
  
  expect_error(link_transform(0.3, list(link = "unknown")))
})

## Edge Cases for link_transform_var() ----
test_that("link_transform_var() handles zero and NA values", {
  expect_error(link_transform_var(0, 100, strategy$family))
  expect_error(link_transform_var(NA, 100, strategy$family))
  expect_error(link_transform_var(-10, 100, strategy$family))
})

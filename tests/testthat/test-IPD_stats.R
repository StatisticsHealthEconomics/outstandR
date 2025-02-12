# from chatgpt

## Setup ----
strategy_maic <- list(
  class = "maic",
  R = 1000,
  formula = y ~ trt,
  family = binomial()
)

strategy_stc <- list(
  class = "stc",
  formula = y ~ trt,
  family = binomial()
)

strategy_gcomp_ml <- list(
  class = "gcomp_ml",
  R = 1000,
  formula = y ~ trt,
  family = binomial()
)

strategy_gcomp_stan <- list(
  class = "gcomp_stan",
  formula = y ~ trt,
  family = binomial()
)

strategy_mim <- list(
  class = "mim",
  formula = y ~ trt,
  family = binomial()
)

ald <- list(
  y.A.sum = 30,
  N.A = 100,
  y.C.sum = 20,
  N.C = 100
)

ipd <- data.frame(
  y = c(1, 0, 1, 0, 1, 0, 1, 0),
  trt = c("A", "A", "A", "A", "C", "C", "C", "C")
)

## General Tests ----
test_that("IPD_stats() works for MAIC", {
  res <- IPD_stats(strategy_maic, ipd, ald)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

test_that("IPD_stats() works for STC", {
  res <- IPD_stats(strategy_stc, ipd, ald)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

test_that("IPD_stats() works for G-computation (ML)", {
  res <- IPD_stats(strategy_gcomp_ml, ipd, ald)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

test_that("IPD_stats() works for G-computation (Stan)", {
  res <- IPD_stats(strategy_gcomp_stan, ipd, ald)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

test_that("IPD_stats() works for Multiple Imputation Marginalisation", {
  res <- IPD_stats(strategy_mim, ipd, ald)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

## Edge Cases ----
test_that("IPD_stats() handles NULL or empty inputs", {
  expect_error(IPD_stats(strategy_maic, NULL, ald))
  expect_error(IPD_stats(strategy_maic, ipd, NULL))
  expect_error(IPD_stats(strategy_maic, list(), ald))
})

test_that("IPD_stats() handles unexpected input types", {
  ipd_wrong <- list(y = "1", trt = "A")
  ald_wrong <- list(y.A.sum = "thirty", N.A = "one hundred")
  
  expect_error(IPD_stats(strategy_maic, ipd_wrong, ald))
  expect_error(IPD_stats(strategy_maic, ipd, ald_wrong))
})

test_that("IPD_stats() handles extreme values", {
  ipd_extreme <- data.frame(
    y = c(1, 1, 1, 1, 0, 0, 0, 0),
    trt = c("A", "A", "A", "A", "C", "C", "C", "C")
  )
  ald_extreme <- list(
    y.A.sum = 0,   # Zero events
    N.A = 100,
    y.C.sum = 100, # All events
    N.C = 100
  )
  res <- IPD_stats(strategy_maic, ipd_extreme, ald_extreme)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

test_that("IPD_stats() handles unsupported strategies", {
  strategy_invalid <- list(class = "unsupported")
  expect_error(IPD_stats(strategy_invalid, ipd, ald))
})

test_that("IPD_stats() handles missing columns", {
  ipd_missing <- data.frame(
    y = c(1, 0, 1, 0),
    # trt column missing
    z = c("A", "A", "C", "C")
  )
  expect_error(IPD_stats(strategy_maic, ipd_missing, ald))
})

test_that("IPD_stats() handles different link functions", {
  strategy_log <- list(class = "stc", formula = y ~ trt, family = binomial(link = "log"))
  strategy_identity <- list(class = "stc", formula = y ~ trt, family = binomial(link = "identity"))
  
  res_log <- IPD_stats(strategy_log, ipd, ald)
  res_identity <- IPD_stats(strategy_identity, ipd, ald)
  
  expect_type(res_log$mean, "double")
  expect_type(res_log$var, "double")
  expect_type(res_identity$mean, "double")
  expect_type(res_identity$var, "double")
})

test_that("IPD_stats() handles unsupported link functions", {
  strategy_unknown <- list(class = "stc", formula = y ~ trt, family = list(link = "unknown"))
  expect_error(IPD_stats(strategy_unknown, ipd, ald))
})

test_that("IPD_stats() handles negative or NA values", {
  ipd_negative <- data.frame(
    y = c(-1, 0, 1, 0),
    trt = c("A", "A", "C", "C")
  )
  ald_na <- list(
    y.A.sum = NA,
    N.A = 100,
    y.C.sum = 20,
    N.C = 100
  )
  expect_error(IPD_stats(strategy_maic, ipd_negative, ald))
  expect_error(IPD_stats(strategy_maic, ipd, ald_na))
})

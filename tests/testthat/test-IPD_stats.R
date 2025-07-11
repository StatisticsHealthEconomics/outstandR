# from chatgpt

library(tibble)

strategy_maic <- list(
  R = 1000,
  formula = y ~ trt,
  trt_var = "trt",
  family = binomial()
) |> 
  `attr<-`(which = "class", value = "maic")

strategy_stc <- list(
  formula = y ~ trt,
  trt_var = "trt",
  family = binomial()
) |> 
  `attr<-`(which = "class", value = "stc")

strategy_gcomp_ml <- list(
  R = 1000,
  formula = y ~ trt,
  trt_var = "trt",
  family = binomial()
) |> 
  `attr<-`(which = "class", value = "gcomp_ml")

strategy_gcomp_stan <- list(
  formula = y ~ trt,
  trt_var = "trt",
  family = binomial()
) |> 
  `attr<-`(which = "class", value = "gcomp_stan")

strategy_mim <- list(
  formula = y ~ trt,
  trt_var = "trt",
  family = binomial()
) |> 
  `attr<-`(which = "class", value = "mim")

ald <- tribble(
  ~variable, ~trt, ~statistic, ~value,
  "y",       "B",  "sum",     30,
  "y",       "C",  "sum",     20,
  NA,        "B",  "N",       100,
  NA,        "C",  "N",       100
)

ipd <- data.frame(
  y = c(1, 0, 1, 0, 1, 0, 1, 0),
  trt = c("A", "A", "A", "A", "C", "C", "C", "C")
)

analysis_params <- 
  list(ipd = ipd,
       ald = ald,
       scale = "log_odds")

## General Tests ----
test_that("calc_IPD_stats() works for MAIC", {
  res <- calc_IPD_stats(strategy = strategy_maic, analysis_params = analysis_params)
  
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

test_that("calc_IPD_stats() works for STC", {
  res <- calc_IPD_stats(strategy_stc, analysis_params)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

test_that("calc_IPD_stats() works for G-computation (ML)", {
  res <- calc_IPD_stats(strategy_gcomp_ml, analysis_params)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

test_that("calc_IPD_stats() works for G-computation (Stan)", {
  res <- calc_IPD_stats(strategy_gcomp_stan, analysis_params)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

test_that("calc_IPD_stats() works for Multiple Imputation Marginalisation", {
  res <- calc_IPD_stats(strategy_mim, analysis_params)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

## Edge Cases ----
# test_that("calc_IPD_stats() handles NULL or empty inputs", {
#   expect_error(calc_IPD_stats(strategy_maic, NULL, ald, scale = "log_odds"))
#   expect_error(calc_IPD_stats(strategy_maic, ipd, NULL, scale = "log_odds"))
#   expect_error(calc_IPD_stats(strategy_maic, list(), ald, scale = "log_odds"))
# })
# 
# test_that("calc_IPD_stats() handles unexpected input types", {
#   ipd_wrong <- list(y = "1", trt = "A")
#   ald_wrong <- list(y.A.sum = "thirty", N.A = "one hundred")
#   
#   expect_error(calc_IPD_stats(strategy_maic, ipd_wrong, ald, scale = "log_odds"))
#   expect_error(calc_IPD_stats(strategy_maic, ipd, ald_wrong, scale = "log_odds"))
# })

test_that("calc_IPD_stats() handles extreme values", {
  ipd_extreme <- data.frame(
    y = c(1, 1, 1, 1, 0, 0, 0, 0),
    trt = c("A", "A", "A", "A", "C", "C", "C", "C")
  )
  
  ald_extreme <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "B",  "sum",     0,     # zero events
    "y",       "C",  "sum",     100,   # all events
    NA,        "B",  "N",       100,
    NA,        "C",  "N",       100
  )
  
  params_extreme <- list(ipd = ipd_extreme,
                         ald = ald_extreme,
                         scale = "log_odds")
  
  res <- calc_IPD_stats(strategy_maic, params_extreme)
  expect_type(res$mean, "double")
  expect_type(res$var, "double")
})

test_that("calc_IPD_stats() handles unsupported strategies", {
  strategy_invalid <- list(class = "unsupported")
  expect_error(calc_IPD_stats(strategy_invalid, analysis_params))
})

test_that("calc_IPD_stats() handles missing columns", {
  ipd_missing <- data.frame(
    y = c(1, 0, 1, 0),
    # trt column missing
    z = c("A", "A", "C", "C")
  )
  
  params_missing <- analysis_params
  params_missing$ipd <- ipd_missing
  
  expect_error(calc_IPD_stats(strategy_maic, params_missing))
})

test_that("calc_IPD_stats() handles unsupported link functions", {
  strategy_unknown <- list(class = "stc",
                           formula = y ~ trt,
                           family = list(link = "unknown"))
  
  expect_error(calc_IPD_stats(strategy_unknown, analysis_params))
})

test_that("calc_IPD_stats() handles negative or NA values", {
  ipd_negative <- data.frame(
    y = c(-1, 0, 1, 0),
    trt = c("A", "A", "C", "C")
  )
  
  ald_na <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "y",       "B",  "sum",     NA,
    "y",       "C",  "sum",     20,
    NA,        "B",  "N",       100,
    NA,        "C",  "N",       100
  )
  
  params_NA_neg <- list(
    ipd = ipd_negative,
    ald = ald_na,
    scale = "log_odds"
  )
  
  expect_error(calc_IPD_stats(strategy_maic, params_NA_neg))
  expect_error(calc_IPD_stats(strategy_maic, params_NA_neg))
})

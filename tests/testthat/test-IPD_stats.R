# IPD_stat tests

library(tibble)

# mock strategy objects ---

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
  N = 1000L,
  formula = y ~ trt,
  trt_var = "trt",
  family = binomial()
) |> 
  `attr<-`(which = "class", value = "gcomp_ml")

strategy_gcomp_bayes <- list(
  formula = y ~ trt,
  trt_var = "trt",
  N = 1000L,
  family = binomial()
) |> 
  `attr<-`(which = "class", value = "gcomp_bayes")

strategy_mim <- list(
  formula = y ~ trt,
  trt_var = "trt",
  N = 1000L,
  family = binomial()
) |> 
  `attr<-`(which = "class", value = "mim")

# aggregate-level outcome data
ald <- tribble(
  ~variable, ~trt, ~statistic, ~value,
  "y",       "B",  "sum",     30,
  "y",       "C",  "sum",     20,
  NA,        "B",  "N",       100,
  NA,        "C",  "N",       100
)

ipd <- data.frame(
  y = sample(c(1, 0), replace = TRUE, size = 40),
  trt = c(rep("A", 20), rep("C", 20))
)

# internal object
analysis_params <- 
  list(ipd = ipd,
       ald = ald,
       ref_trt = "C",
       ipd_comp = "A",
       scale = "log_odds")

## test when no covariates

test_that("calc_IPD_stats() works for MAIC", {

  res <- calc_IPD_stats(strategy_maic, analysis_params)
  
  expect_type(res$contrasts$mean, "double")
  expect_type(res$contrasts$var, "double")

  # single arm ipd
  
  params_single_ipd <- list(
    ipd = 
      data.frame(variable = "y",
                 trt = "B",
                 statistic = "sum",
                 value = 30),
    ald = ald,
    scale = "log_odds")
  
  strategy_maic_single_sample <- list(
    R = 1,   # results in TWO samples, including original
    formula = y ~ trt,
    trt_var = "trt",
    family = binomial()
  ) |> 
    `attr<-`(which = "class", value = "maic")
  
  calc_IPD_stats(strategy_maic_single_sample, params_single_ipd) |> 
    expect_warning(regexp = "Bootstrap sample contains less than two treatment levels. Returning NA.") |>  
    expect_warning(regexp = "Bootstrap sample contains less than two treatment levels. Returning NA.") 
})

test_that("calc_IPD_stats() works for STC", {
  res <- calc_IPD_stats(strategy_stc, analysis_params)
  
  expect_type(res$contrasts$mean, "double")
  expect_type(res$contrasts$var, "double")
})

test_that("calc_IPD_stats() works for G-computation (ML)", {
  res_gcomp_ml_null <- calc_IPD_stats(strategy_gcomp_ml, analysis_params)
  expect_length(res_gcomp_ml_null, 4)
})

test_that("calc_IPD_stats() works for G-computation (Stan)", {
    res_gcomp_bayes_null <- calc_IPD_stats(strategy_gcomp_bayes, analysis_params)
    expect_length(res_gcomp_bayes_null, 4)
})

test_that("calc_IPD_stats() works for Multiple Imputation Marginalisation", {
  res_mim_null <- calc_IPD_stats(strategy_mim, analysis_params)
  expect_length(res_mim_null, 4)
})

## edge cases

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
  
  res <- suppressWarnings(
    calc_IPD_stats(strategy_maic, params_extreme)
  )
  
  expect_type(res$contrasts$mean, "double")
  expect_type(res$contrasts$var, "double")
})

test_that("calc_IPD_stats() handles unsupported strategies", {
  strategy_invalid <- list() |> 
    `attr<-`(which = "class",
             value = "unsupported")
  
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
  
  # calc_IPD_stats(strategy_stc, params_missing)
  # 
  # expect_error()
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

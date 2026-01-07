# test strategies maic, stc

library(dplyr)
library(glue)


test_that("different combinations of covariates in formula", {
  
  load(test_path("testdata/BC_ALD.RData"))
  load(test_path("testdata/AC_IPD.RData"))
  
  # maic
  expect_error(strategy_maic(formula = as.formula("y ~ 1")),
               regexp = "Treatment term 'trt' is missing in the formula")

  expect_message(strategy_maic(formula = as.formula("y ~ X3 + X4")),
                 regexp = "Treatment is guessed as:")
  
  strat_1234 <- strategy_maic(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  strat_123 <- strategy_maic(formula = as.formula("y ~ X1 + X2 + trt + trt:X2 + trt:X3"))
  strat_31 <- strategy_maic(formula = as.formula("y ~ X3 + trt*X1"))
  strat_13 <- strategy_maic(formula = as.formula("y ~ trt*X1 + X3"))
  strat_1 <- strategy_maic(formula = as.formula("y ~ trt*X1"))
  
  res <- outstandR(AC_IPD, BC_ALD, strategy = strat_1234)
  expect_length(res, 11)
  
  res <- outstandR(AC_IPD, BC_ALD, strategy = strat_123)
  
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_31))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_13))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_1))

  # stc
  expect_error(strategy_stc(formula = as.formula("y ~ 1")),
               regexp = "Treatment term 'trt' is missing in the formula")
  
  expect_message(strategy_stc(formula = as.formula("y ~ X3 + X4")),
                 regexp = "Treatment is guessed as:")
  
  strat_1234 <- strategy_stc(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  strat_31 <- strategy_stc(formula = as.formula("y ~ X3 + trt*X1"))
  strat_13 <- strategy_stc(formula = as.formula("y ~ trt*X1 + X3"))
  strat_1 <- strategy_stc(formula = as.formula("y ~ trt*X1"))

  expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_1234)$results$contrasts$means$AC,
               expected = -0.27, tolerance = 0.1)
  
  expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_31)$results$contrasts$means$AC,
               expected = -0.27, tolerance = 0.1)
  
  expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_13)$results$contrasts$means$AC,
               expected = -0.27, tolerance = 0.1)
  
  expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_1)$results$contrasts$means$AC,
               expected = -0.30, tolerance = 0.1)
})

test_that("compare with maicplus package with binary outcome", {
  
  ## original
  
  library(maicplus)
  
  data(centered_ipd_twt, package = "maicplus")
  data(adrs_twt, package = "maicplus")
  data(adsl_sat, package = "maicplus")  # original IPD
  
  adsl_sat <- adsl_sat |> 
    mutate(SEX_MALE = as.integer(SEX == "Male"))
  
  # single continuous and binary covariate case
  adsl_colnames <- c("AGE", "SEX_MALE")
  centered_colnames <- paste0(adsl_colnames, "_CENTERED")
  
  weighted_data <- estimate_weights(
    data = centered_ipd_twt,
    centered_colnames = centered_colnames)
  
  agd <- data.frame(
    AGE_MEAN = 51,
    SEX_MALE_PROP = 147 / 300)
  
  binary_agd <- data.frame(
    ARM = rep(c("B", "C"), each = 2),
    RESPONSE = rep(c("YES", "NO"), 2),
    COUNT = c(280, 200, 120, 200))
  
  pseudo_adrs <- get_pseudo_ipd_binary(
    binary_agd = binary_agd,
    format = "stacked")
  
  res_maicplus <-
    maic_anchored(
      weights_object = weighted_data,
      ipd = adrs_twt,
      pseudo_ipd = pseudo_adrs,
      trt_ipd = "A",
      trt_agd = "B",
      trt_common = "C",
      normalize_weight = FALSE,
      endpoint_type = "binary",
      endpoint_name = "Binary Endpoint",
      eff_measure = "OR")
  
  ## {outstandR}
  
  lin_form <- as.formula(glue::glue("y ~ trt * ({paste(adsl_colnames, collapse = ' + ')})"))
  
  AC.IPD <- adsl_twt |>
    merge(adrs_twt) |> 
    rename(trt = ARM,
           y = RESPONSE)
  
  BC.ALD <- tribble(
    ~variable, ~trt, ~statistic, ~value,
    "AGE",      NA,   "mean",     agd$AGE_MEAN,
    "SEX_MALE", NA,   "prop",     agd$SEX_MALE_PROP,
    "y",        "B",  "sum",      binary_agd[binary_agd$ARM == "B" &
                                             binary_agd$RESPONSE == "YES", "COUNT"],
    "y",        "C",  "sum",      binary_agd[binary_agd$ARM == "C" &
                                             binary_agd$RESPONSE == "YES", "COUNT"],
    NA,         "B",  "N",        480,
    NA,         "C",  "N",        320
  )
  
  res_outstandr <- 
    maic.boot(ipd = AC.IPD,
              formula = lin_form,
              family = binomial("logit"),
              ald = BC.ALD,
              trt_var = "trt")

  res_outstandr_unadjusted <- 
    maic.boot(ipd = AC.IPD,
              formula = lin_form,
              family = binomial("logit"),
              ald = BC.ALD,
              trt_var = "trt",
              hat_w = rep(1, nrow(AC.IPD)))
  
  maicplus_AC <-
    with(res_maicplus$inferential,
         summary[summary$case == "adjusted_AC", "OR"])

  maicplus_AC_unadjusted <-
    with(res_maicplus$inferential,
         summary[summary$case == "AC", "OR"])
  
  outstandr_OR <- with(as.list(res_outstandr), pC/(1-pC)/(pA/(1-pA)))
  outstandr_OR_unadjusted <- with(as.list(res_outstandr_unadjusted), pC/(1-pC)/(pA/(1-pA)))
  
  expect_equal(maicplus_AC, outstandr_OR)
  expect_equal(maicplus_AC_unadjusted, outstandr_OR_unadjusted)
  
  
  # different output scales
  ##TODO:
  # calculate_ate()
  
})

test_that("compare with maicplus package with continuous outcome", {

  n_ipd <- 100
  
  ipd <- data.frame(
    age = rnorm(n_ipd, mean = 55, sd = 10),
    sex = rbinom(n_ipd, 1, 0.5),  # 1 = male
    y = rnorm(n_ipd, mean = 60, sd = 8)
  )
  
  agd_covariates <- list(age = 60, sex = 0.6)
  
  EM_mat <- select(ipd, - y)
  EM_mat$age <- EM_mat$age - agd_covariates$age
  EM_mat$sex <- EM_mat$sex - agd_covariates$sex
  
  weights <- maic_weights(X_EM = EM_mat)
  
  ipd$weight <- weights
  weighted_mean_ipd <- weighted.mean(ipd$y, ipd$weight)
})


#
test_that("mismatch between covariates in ald and ipd / formula", {
  
  # different order between ald and ipd
  # age sex; sex age
  
  # different order between ald and formula
  # age + sex; mean.sex mean.age
  
  # same covariates as EM and PF
  # age + trt * (age)
  
  # transformed covariates in formula
  # I(X1^2)
})


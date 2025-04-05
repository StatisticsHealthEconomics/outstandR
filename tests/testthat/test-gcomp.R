#

library(dplyr)

#
test_that("different combinations of covariates in formula", {
  
  load(test_path("testdata/BC_ALD.RData"))
  load(test_path("testdata/AC_IPD.RData"))
  
  # gcomp_ml
  
  expect_error(strategy_gcomp_ml(formula = as.formula("y ~ 1")),
               regexp = "Treatment term, trt, is missing in the formula")
  
  expect_error(strategy_gcomp_ml(formula = as.formula("y ~ X3 + X4")),
               regexp = "Treatment term, trt, is missing in the formula")
  
  strat_1234 <- strategy_gcomp_ml(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  strat_31 <- strategy_gcomp_ml(formula = as.formula("y ~ X3 + trt*X1"))
  strat_13 <- strategy_gcomp_ml(formula = as.formula("y ~ trt*X1 + X3"))
  strat_1 <- strategy_gcomp_ml(formula = as.formula("y ~ trt*X1"))
  
  expect_length(outstandR(AC_IPD, BC_ALD, strategy = strat_1234), 3)
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_31))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_13))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_1))
  
  # gcomp_stan
  
  expect_error(strategy_gcomp_stan(formula = as.formula("y ~ 1")),
               regexp = "Treatment term, trt, is missing in the formula")
  
  expect_error(strategy_gcomp_stan(formula = as.formula("y ~ X3 + X4")),
               regexp = "Treatment term, trt, is missing in the formula")
  
  strat_1234 <- strategy_gcomp_stan(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  strat_31 <- strategy_gcomp_stan(formula = as.formula("y ~ X3 + trt*X1"))
  strat_13 <- strategy_gcomp_stan(formula = as.formula("y ~ trt*X1 + X3"))
  strat_1 <- strategy_gcomp_stan(formula = as.formula("y ~ trt*X1"))
  
  expect_length(outstandR(AC_IPD, BC_ALD, strategy = strat_1234), 3)
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_31))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_13))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_1))
})


test_that("compare with stdReg2 package for continuous outcome", {
  
  # original
  library(stdReg2)
  library(causaldata)
  
  # load dataset
  nhefs_dat <- causaldata::nhefs_complete
  
  # m <- glm(wt82_71 ~ qsmk + sex + age, 
  #          data = nhefs_dat)
  
  res_stdReg2 <-
    standardize_glm(
      wt82_71 ~ qsmk + sex + age,
      data = nhefs_dat, 
      values = list(qsmk = c(0,1)),
      contrasts = c("difference", "ratio"),
      reference = 0)
  
  ## {outstandR}
  
  nhefs_ipd <- nhefs_dat |> 
    select(qsmk, sex, age, wt82_71) |> 
    rename(trt = qsmk,
           y = wt82_71) |> 
    mutate(sex = as.numeric(sex) - 1)
  
  lin_form <- as.formula(y ~ trt * (sex + age))
  
  # create aggregate data
  nhefs.X <- nhefs_ipd |> 
    summarise(across(-c(trt, y),
                     list(mean = mean, sd = sd),
                     .names = "{fn}.{col}"))
  
  nhefs.B <- dplyr::filter(nhefs_ipd, trt == 1) |> 
    summarise(y.B.sum = sum(y),
              y.B.bar = mean(y),
              y.B.sd = sd(y),
              N.B = n())
  
  nhefs.C <- dplyr::filter(nhefs_ipd, trt == 0) |> 
    summarise(y.C.sum = sum(y),
              y.C.bar = mean(y),
              y.C.sd = sd(y),
              N.C = n())
  
  nhefs_ald <- cbind.data.frame(nhefs.X, nhefs.C, nhefs.B)
  
  res_outstandr <- gcomp_ml_means(
    formula = lin_form,
    family = "gaussian",
    ald = nhefs_ald,
    ipd = nhefs_ipd,
    N = 100000)
  
  # means
  
  stdReg2_estimates <- res_stdReg2$res$estimates
  
  expect_equal(unname(res_outstandr["0"]),
               stdReg2_estimates$estimates[stdReg2_estimates$qsmk == 0],
               tolerance = 0.1)
  
  expect_equal(unname(res_outstandr["1"]),
               stdReg2_estimates$estimates[stdReg2_estimates$qsmk == 1],
               tolerance = 0.1)
  
  # std errors
  ##TODO:
  
  # different output scales
  
  expect_equal(
    res_stdReg2$res_contrast[[2]]$est_table$Estimate[2],
    calculate_ate(mean_A = unname(res_outstandr["1"]),
                  mean_C = unname(res_outstandr["0"]),
                  effect = "risk_difference"),
    tolerance = 0.01)
  
  expect_equal(
    res_stdReg2$res_contrast[[3]]$est_table$Estimate[2],
    exp(calculate_ate(mean_A = unname(res_outstandr["1"]),
                      mean_C = unname(res_outstandr["0"]),
                      effect = "log_relative_risk")),
    tolerance = 0.01)
  
})

test_that("compare with stdReg2 package for binary outcome", {
  
  # Step 1: Simulate data
  n <- 1000
  data <- data.frame(
    age = rnorm(n, mean = 50, sd = 10),
    sex = rbinom(n, 1, 0.5),
    trt = rbinom(n, 1, 0.5)
  )
  
  # Generate a binary outcome with logit link
  logit_p <- -2 + 0.05 * data$age + 0.5 * data$sex + 1.2 * data$treat
  prob <- 1 / (1 + exp(-logit_p))
  data$y <- rbinom(n, 1, prob)
  
  # Step 2: Fit logistic regression
  fit <- glm(y ~ trt + age + sex, data = data, family = binomial)
  
  # Step 3: Standardize to whole sample (marginal treatment effect)
  res_stdReg2 <- stdGlm(
    fit = fit,
    data = data,
    X = "trt",
    x = c(0, 1),  # treatment levels to compare
    outcome.type = "binary"
  )
  
  ## {outstandR}
  
  # create aggregate data
  nhefs.X <- data |> 
    summarise(across(-c(trt, y),
                     list(mean = mean, sd = sd),
                     .names = "{fn}.{col}"))
  
  nhefs.B <- dplyr::filter(data, trt == 1) |> 
    summarise(y.B.sum = sum(y),
              y.B.bar = mean(y),
              y.B.sd = sd(y),
              N.B = n())
  
  nhefs.C <- dplyr::filter(data, trt == 0) |> 
    summarise(y.C.sum = sum(y),
              y.C.bar = mean(y),
              y.C.sd = sd(y),
              N.C = n())
  
  data_ald <- cbind.data.frame(nhefs.X, nhefs.C, nhefs.B)
  
  lin_form <- as.formula(y ~ trt * (sex + age))
  
  res_outstandr <- gcomp_ml_means(formula = lin_form,
                                  family = "binomial",
                                  ald = data_ald,
                                  ipd = data,
                                  N = 10000)
  # means
  
  stdReg2_estimates <- res_stdReg2$res$estimates
  
  expect_equal(unname(res_outstandr["0"]),
               stdReg2_estimates$estimates[stdReg2_estimates$qsmk == 0],
               tolerance = 0.1)
  
  expect_equal(unname(res_outstandr["1"]),
               stdReg2_estimates$estimates[stdReg2_estimates$qsmk == 1],
               tolerance = 0.1)
  
  # std errors
  ##TODO:
  
  # different output scales
  
  expect_equal(
    res_stdReg2$res_contrast[[2]]$est_table$Estimate[2],
    calculate_ate(mean_A = unname(res_outstandr["1"]),
                  mean_C = unname(res_outstandr["0"]),
                  effect = "risk_difference"),
    tolerance = 0.01)
  
  expect_equal(
    res_stdReg2$res_contrast[[3]]$est_table$Estimate[2],
    exp(calculate_ate(mean_A = unname(res_outstandr["1"]),
                      mean_C = unname(res_outstandr["0"]),
                      effect = "log_relative_risk")),
    tolerance = 0.01)
  
})

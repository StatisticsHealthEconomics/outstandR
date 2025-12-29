# G-computation

library(dplyr)
library(stdReg2)

#
test_that("different combinations of covariates in formula", {
  
  load(test_path("testdata/BC_ALD.RData"))
  load(test_path("testdata/AC_IPD.RData"))
  
  BC_ALD <- reshape_ald_to_long(BC_ALD)
  AC_IPD$trt <- factor(AC_IPD$trt, labels = c("C", "A"))  # from 0, 1
  
  ### gcomp_ml
  
  expect_error(strategy_gcomp_ml(formula = as.formula("y ~ 1")),
               regexp = "Treatment term 'trt' is missing in the formula")
  
  expect_error(strategy_gcomp_ml(formula = as.formula("y ~ X3 + X4"), trt_var = "trt"),
               regexp = "Treatment term 'trt' is missing in the formula")
  
  strat_1234 <- strategy_gcomp_ml(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  strat_31 <- strategy_gcomp_ml(formula = as.formula("y ~ X3 + trt*X1"))
  strat_13 <- strategy_gcomp_ml(formula = as.formula("y ~ trt*X1 + X3"))
  strat_1 <- strategy_gcomp_ml(formula = as.formula("y ~ trt*X1"))
  
  res <- outstandR(ipd_trial = AC_IPD, ald_trial = BC_ALD, strategy = strat_1234)
  
  expect_length(res, 2)
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_31))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_13))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_1))
  
  ### gcomp_bayes
  
  expect_error(strategy_gcomp_bayes(formula = as.formula("y ~ 1")),
               regexp = "Treatment term 'trt' is missing in the formula")
  
  expect_error(strategy_gcomp_bayes(formula = as.formula("y ~ X3 + X4"), trt_var = "trt"),
               regexp = "Treatment term 'trt' is missing in the formula")
  
  strat_1234 <- strategy_gcomp_bayes(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"))
  strat_31 <- strategy_gcomp_bayes(formula = as.formula("y ~ X3 + trt*X1"))
  strat_13 <- strategy_gcomp_bayes(formula = as.formula("y ~ trt*X1 + X3"))
  strat_1 <- strategy_gcomp_bayes(formula = as.formula("y ~ trt*X1"))
  
  res <- outstandR(ipd_trial = AC_IPD, ald_trial = BC_ALD, strategy = strat_1234)
  
  expect_length(res, 2)
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_31))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_13))
  # expect_equal(outstandR(AC_IPD, BC_ALD, strategy = strat_1))
})


test_that("compare with stdReg2 package for continuous outcome", {
  
  ## original
  
  library(causaldata)
  
  # load dataset
  nhefs_dat <- causaldata::nhefs_complete
  
  # m <- glm(wt82_71 ~ qsmk + sex + age, 
  #          data = nhefs_dat)
  
  res_stdReg2 <-
    stdReg2::standardize_glm(
      wt82_71 ~ qsmk + sex + age,
      data = nhefs_dat, 
      values = list(qsmk = c(0,1)),
      family = "gaussian",
      contrasts = c("difference", "ratio"),
      reference = 0)
  
  ## {outstandR}
  
  nhefs_ipd <- nhefs_dat |> 
    dplyr::select(qsmk, sex, age, wt82_71) |> 
    dplyr::rename(trt = qsmk,
                  y = wt82_71) |> 
    dplyr::mutate(sex = as.numeric(sex) - 1,
                  trt = factor(trt, labels = c("C", "A")))
  
  lin_form <- as.formula(y ~ trt * (sex + age))
  
  # create aggregate data
  nhefs.X <- nhefs_ipd |> 
    dplyr::summarise(across(-c(trt, y),
                            list(mean = mean, sd = sd),
                            .names = "{fn}.{col}")) |> 
    reshape2::melt() |> 
    tidyr::separate(variable,
                    into = c("statistic", "variable"),
                    sep = "\\.") |> 
    dplyr::mutate(trt = NA)
  
  ald_out <- nhefs_ipd |> 
    group_by(trt) |> 
    dplyr::summarise(
      y.sum = sum(y),
      y.mean = mean(y),
      y.sd = sd(y),
      .N = dplyr::n()) |> 
    reshape2::melt(id.vars = "trt") |> 
    tidyr::separate(variable,
                    into = c("variable", "statistic"),
                    sep = "\\.") |> 
    dplyr::mutate(trt = ifelse(trt == "A", "B", "C"))
  
  nhefs_ald <- 
    dplyr::bind_rows(nhefs.X, ald_out) |>
    tibble::as_tibble() |> 
    dplyr::mutate(trt = factor(trt))
  
  res_outstandr <- gcomp_ml_means(
    formula = lin_form,
    family = "gaussian",
    ipd = nhefs_ipd,
    ald = nhefs_ald,
    ref_trt = "C",
    comp_trt = "A",
    trt_var = "trt",
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
    calculate_ate(mean_comp = unname(res_outstandr["1"]),
                  mean_ref = unname(res_outstandr["0"]),
                  effect = "risk_difference"),
    tolerance = 0.01)
  
  expect_equal(
    res_stdReg2$res_contrast[[3]]$est_table$Estimate[2],
    exp(calculate_ate(mean_comp = unname(res_outstandr["1"]),
                      mean_ref = unname(res_outstandr["0"]),
                      effect = "log_relative_risk")),
    tolerance = 0.01)
})

test_that("compare with stdReg2 package for binary outcome", {
  
  # simulate data
  n <- 1000
  
  data <- data.frame(
    age = rnorm(n, mean = 50, sd = 10),
    sex = rbinom(n, 1, 0.5),
    trt = rbinom(n, 1, 0.5))
  
  # generate a binary outcome with logit link
  logit_p <- -1 + data$trt * (0.03 * data$age + 0.2 * data$sex)
  prob <- 1 / (1 + exp(-logit_p))
  data$y <- rbinom(n, 1, prob)
  
  res_stdReg2 <-
    standardize_glm(
      y ~ trt + age + sex,
      family = "binomial",
      data = data, 
      values = list(trt = c(0,1)),
      contrasts = c("difference", "ratio"),
      reference = 0)
  
  ## {outstandR}
  
  # create aggregate data
  data.X <- data |> 
    dplyr::summarise(across(-c(trt, y),
                            list(mean = mean, sd = sd),
                            .names = "{fn}.{col}")) |> 
    reshape2::melt() |> 
    tidyr::separate(variable,
                    into = c("statistic", "variable"),
                    sep = "\\.") |> 
    dplyr::mutate(trt = NA)
  
  ald_out <- data |> 
    group_by(trt) |> 
    dplyr::summarise(
      y.sum = sum(y),
      y.mean = mean(y),
      y.sd = sd(y),
      .N = dplyr::n()) |> 
    reshape2::melt(id.vars = "trt") |> 
    tidyr::separate(variable,
                    into = c("variable", "statistic"),
                    sep = "\\.") |> 
    dplyr::mutate(trt = ifelse(trt == "1", "B", "C"))
  
  data_ald <- bind_rows(data.X, ald_out) |>
    tibble::as_tibble() |> 
    dplyr::mutate(trt = factor(trt))
  
  data$trt <- factor(data$trt, labels = c("C", "A"))
  
  lin_form <- as.formula(y ~ trt * (sex + age))
  
  res_outstandr <- gcomp_ml_means(
    formula = lin_form,
    family = "binomial",
    ipd = data,
    ald = data_ald,
    trt_var = "trt",
    comp_trt = "A",
    ref_trt = "C",
    N = 10000)
  
  # means
  
  stdReg2_estimates <- res_stdReg2$res$estimates
  
  expect_equal(unname(res_outstandr["0"]),
               stdReg2_estimates$estimates[stdReg2_estimates$trt == 0],
               tolerance = 0.1)
  
  expect_equal(unname(res_outstandr["1"]),
               stdReg2_estimates$estimates[stdReg2_estimates$trt == 1],
               tolerance = 0.1)
  
  # std errors
  ##TODO:
  
  # different output scales
  
  expect_equal(
    res_stdReg2$res_contrast[[2]]$est_table$Estimate[2],
    calculate_ate(mean_comp = unname(res_outstandr["1"]),
                  mean_ref = unname(res_outstandr["0"]),
                  effect = "risk_difference"),
    tolerance = 0.01)
  
  expect_equal(
    res_stdReg2$res_contrast[[3]]$est_table$Estimate[2],
    exp(calculate_ate(mean_comp = unname(res_outstandr["1"]),
                      mean_ref = unname(res_outstandr["0"]),
                      effect = "log_relative_risk")),
    tolerance = 0.01)
})


test_that("compare with marginaleffects package for binary outcome", {
  
  # simulate data
  n <- 1000
  
  data <- data.frame(
    age = rnorm(n, mean = 50, sd = 10),
    sex = rbinom(n, 1, 0.5),
    trt = rbinom(n, 1, 0.5))
  
  # generate a binary outcome with logit link
  logit_p <- -1 + data$trt * (0.03 * data$age + 0.2 * data$sex)
  prob <- 1 / (1 + exp(-logit_p))
  data$y <- rbinom(n, 1, prob)
  data$trt <- factor(data$trt, labels = c("C", "A"))
  
  lin_form <- as.formula(y ~ trt * (sex + age))
  
  fit <- glm(lin_form, family = binomial(link = "logit"), data = data)
  
  m_eff_pred <- marginaleffects::avg_predictions(fit, variables = "trt", by = "trt")
  
  ## {outstandR}
  
  # create aggregate data
  data.X <- data |> 
    dplyr::summarise(across(-c(trt, y),
                            list(mean = mean, sd = sd),
                            .names = "{fn}.{col}")) |> 
    reshape2::melt() |> 
    tidyr::separate(variable,
                    into = c("statistic", "variable"),
                    sep = "\\.") |> 
    dplyr::mutate(trt = NA)
  
  ald_out <- data |> 
    group_by(trt) |> 
    dplyr::summarise(
      y.sum = sum(y),
      y.mean = mean(y),
      y.sd = sd(y),
      .N = dplyr::n()) |> 
    reshape2::melt(id.vars = "trt") |> 
    tidyr::separate(variable,
                    into = c("variable", "statistic"),
                    sep = "\\.") |> 
    dplyr::mutate(trt = ifelse(trt == "1", "B", "C"))
  
  data_ald <- bind_rows(data.X, ald_out) |>
    tibble::as_tibble() |> 
    dplyr::mutate(trt = factor(trt))
  
  data$trt <- factor(data$trt, labels = c("C", "A"))
  
  res_outstandr <- gcomp_ml_means(
    formula = lin_form,
    family = "binomial",
    ipd = data,
    ald = data_ald,
    trt_var = "trt",
    comp_trt = "A",
    ref_trt = "C",
    N = 10000)
  
  # means
  
  estimates <- m_eff_pred$estimate
  
  expect_equal(unname(res_outstandr["0"]),
               estimates[1],
               tolerance = 0.1)
  
  expect_equal(unname(res_outstandr["1"]),
               estimates[2],
               tolerance = 0.1)
})


test_that("compare with marginaleffects package for continuous outcome", {
  ## original
  
  library(causaldata)
  
  # load dataset
  nhefs_dat <- causaldata::nhefs_complete
  
  fit <- glm(wt82_71 ~ qsmk * (sex + age), family = "gaussian", data = nhefs_dat)
  
  m_eff_pred <- marginaleffects::avg_predictions(fit, variables = "qsmk", by = "qsmk")
  
  ## {outstandR}
  
  nhefs_ipd <- nhefs_dat |> 
    select(qsmk, sex, age, wt82_71) |> 
    rename(trt = qsmk,
           y = wt82_71) |> 
    mutate(sex = as.numeric(sex) - 1,
           trt = factor(trt, labels = c("C", "A")))
  
  lin_form <- as.formula(y ~ trt * (sex + age))
  
  # create aggregate data
  
  nhefs.X <- nhefs_ipd |> 
    dplyr::summarise(across(-c(trt, y),
                            list(mean = mean, sd = sd),
                            .names = "{fn}.{col}")) |> 
    reshape2::melt() |> 
    tidyr::separate(variable,
                    into = c("statistic", "variable"),
                    sep = "\\.") |> 
    dplyr::mutate(trt = NA)
  
  ald_out <- nhefs_ipd |> 
    group_by(trt) |> 
    dplyr::summarise(
      y.sum = sum(y),
      y.mean = mean(y),
      y.sd = sd(y),
      .N = dplyr::n()) |> 
    reshape2::melt(id.vars = "trt") |> 
    tidyr::separate(variable,
                    into = c("variable", "statistic"),
                    sep = "\\.") |> 
    dplyr::mutate(trt = ifelse(trt == "1", "B", "C"))
  
  nhefs_ald <- bind_rows(nhefs.X, ald_out) |>
    tibble::as_tibble() |> 
    dplyr::mutate(trt = factor(trt))
  
  res_outstandr <- gcomp_ml_means(
    formula = lin_form,
    family = "gaussian",
    ipd = nhefs_ipd,
    ald = nhefs_ald,
    trt_var = "trt",
    comp_trt = "A",
    ref_trt = "C",
    N = 100000)
  
  # means
  
  expect_equal(unname(res_outstandr["0"]),
               m_eff_pred$estimate[1],
               tolerance = 0.1)
  
  expect_equal(unname(res_outstandr["1"]),
               m_eff_pred$estimate[2],
               tolerance = 0.1)
  
  # std errors
  ##TODO:
  
  # different output scales
  
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

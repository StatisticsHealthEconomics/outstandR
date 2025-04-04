
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


test_that("compare with stdReg2 package", {
  
  # original
  library(stdReg2)
  library(causaldata)
  
  # load dataset
  nhefs_dat <- causaldata::nhefs_complete
  
  m <- glm(wt82_71 ~ qsmk + sex + race + age + I(age^2) + as.factor(education) +
             smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
             as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2), 
           data = nhefs_dat)
  
  res_stdReg2 <-
    standardize_glm(
      wt82_71 ~ qsmk + sex + race + age + I(age^2) + as.factor(education) +
        smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
        as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2),
      data = nhefs_dat, 
      values = list(qsmk = c(0,1)),
      contrasts = c("difference", "ratio"),
      reference = 0)
  
  ## {outstandR}
  
  nhefs_ipd <- nhefs_dat |> 
    select(qsmk, sex, race, age, smokeintensity, smokeyrs, wt71,
           wt82_71, education, exercise, active) |> 
    rename(trt = qsmk,
           y = wt82_71) |> 
    mutate(sex = as.numeric(sex) - 1,
           race = as.numeric(race) - 1) |> 
    maicplus::dummize_ipd(dummize_cols = c("education", "exercise", "active"),
                          dummize_ref_level = c("1","0","0")) |> 
    select(-education, -exercise, -active)
  
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
  
  lin_form <-
    as.formula(y ~ trt * (sex + race + age + I(age^2) +
                            smokeintensity + I(smokeintensity^2) +
                            EDUCATION_2 + EDUCATION_3 + EDUCATION_4 + EDUCATION_5 +
                            EXERCISE_1 + EXERCISE_2 + ACTIVE_1 + ACTIVE_2 +
                            smokeyrs + I(smokeyrs^2) + wt71 + I(wt71^2)))
  nhefs_strat <- 
    strategy_gcomp_ml(formula = lin_form,
                      family = gaussian(link = "identity"))
  
  res_outstandr <- outstandR(BC.ALD = nhefs_ald,
                             AC.IPD = nhefs_ipd,
                             strategy = nhefs_strat,
                             scale = "risk_difference")
  
  expect_equal(res_outstandr, res_stdReg2)
})

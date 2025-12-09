#

test_that("simulate_ALD_pseudo_pop directly", {
  
  form <- y ~ PF_cont_1 + PF_cont_2 + trt*EM_cont_1 + trt*EM_cont_2
  corX <- 0.2
  
  ####################
  # without marginals
  
  # providing competing arguments
  
  # no rho ie taken from ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial_cont,
    ald = ald_trial_cont,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "EM_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "EM_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "PF_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "PF_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$EM_cont_1),
    expected = EM1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$EM_cont_2),
    expected = EM2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_cont_1),
    expected = PF1_mean,
    tolerance = 0.01
  )
  expect_equal(
    object = mean(result$PF_cont_2),
    expected = PF2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$EM_cont_1, result$EM_cont_2),
    expected = cor(ipd_trial_cont$EM_cont_1, ipd_trial_cont$EM_cont_2),
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$PF_cont_1, result$PF_cont_2),
    expected = cor(ipd_trial_cont$PF_cont_1, ipd_trial_cont$PF_cont_2),
    tolerance = 0.1
  )
  
  # rho (same) and ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial_cont,
    ald = ald_trial_cont,
    rho = corX,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial_cont |>
    dplyr::filter(variable == "EM_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "EM_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "PF_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "PF_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$EM_cont_1),
    expected = EM1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$EM_cont_2),
    expected = EM2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_cont_1),
    expected = PF1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_cont_2),
    expected = PF2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$PF_cont_1, result$PF_cont_2),
    expected = corX,
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$EM_cont_1, result$EM_cont_2),
    expected = corX,
    tolerance = 0.1
  )
  
  # rho (different) and ipd
  
  rho_new <- 0.6
  
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial_cont,
    ald = ald_trial_cont,
    rho = rho_new,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "EM_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial_cont |>
    dplyr::filter(variable == "EM_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial_cont |>
    dplyr::filter(variable == "PF_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "PF_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$EM_cont_1),
    expected = EM1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$EM_cont_2),
    expected = EM2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_cont_1),
    expected = PF1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_cont_2),
    expected = PF2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$EM_cont_1, result$EM_cont_2),
    expected = rho_new,
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$PF_cont_1, result$PF_cont_2),
    expected = rho_new,
    tolerance = 0.1
  )
  
  # no ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ald = ald_trial_cont,
    rho = rho_new,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "EM_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial_cont |>
    dplyr::filter(variable == "EM_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "PF_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial_cont |> 
    dplyr::filter(variable == "PF_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$EM_cont_1),
    expected = EM1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$EM_cont_2),
    expected = EM2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_cont_1),
    expected = PF1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_cont_2),
    expected = PF2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$EM_cont_1, result$EM_cont_2),
    expected = rho_new,
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$PF_cont_1, result$PF_cont_2),
    expected = rho_new,
    tolerance = 0.1
  )
  
  #################
  # with marginals
  
  meanX_BC <- c(0.6, 0.6)  # mean of each normally-distributed covariate in BC
  sdX <- c(0.4, 0.4)  # standard deviation of each covariate (same for AC and BC)
  meanX_EM_BC <- c(0.6, 0.6)  # mean of each normally-distributed EM covariate in BC
  sdX_EM <- c(0.4, 0.4)  # standard deviation of each EM covariate
  
  marginals_orig <- list(
    marginal_dists = c(EM_cont_1 = "norm",
                       EM_cont_2 = "norm",
                       PF_cont_1 = "norm",
                       PF_cont_2 = "norm"),
    marginal_params = list(
      EM_cont_1 = list(mean = meanX_EM_BC[1], sd = sdX_EM[1]),
      EM_cont_2 = list(mean = meanX_EM_BC[2], sd = sdX_EM[2]),
      PF_cont_1 = list(mean = meanX_BC[1], sd = sdX[1]),
      PF_cont_2 = list(mean = meanX_BC[2], sd = sdX[2])
    )
  )
  
  marginals_new <- list(
    marginal_dists = c(EM_cont_1 = "norm",
                       EM_cont_2 = "norm",
                       PF_cont_1 = "norm",
                       PF_cont_2 = "norm"),
    marginal_params = list(
      EM_cont_1 = list(mean = 0, sd = 1),
      EM_cont_2 = list(mean = 0, sd = 1),
      PF_cont_1 = list(mean = 0, sd = 1),
      PF_cont_2 = list(mean = 0, sd = 1)
    )
  )
  
  # no rho ie taken from ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial_cont,
    ald = ald_trial_cont,
    trt_var = "trt",
    N = 100000,
    marginal_distns = marginals_orig$marginal_dists,
    marginal_params = marginals_orig$marginal_params
  )
  
  expect_equal(
    object = mean(result$EM_cont_1),
    expected = marginals_orig$marginal_params$EM_cont_1$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$EM_cont_2),
    expected = marginals_orig$marginal_params$EM_cont_2$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_cont_1),
    expected = marginals_orig$marginal_params$PF_cont_1$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_cont_2),
    expected = marginals_orig$marginal_params$PF_cont_2$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$EM_cont_1, result$EM_cont_2),
    expected = cor(ipd_trial_cont$EM_cont_1, ipd_trial_cont$EM_cont_2),
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$PF_cont_1, result$PF_cont_2),
    expected = cor(ipd_trial_cont$PF_cont_1, ipd_trial_cont$PF_cont_2),
    tolerance = 0.1
  )
  
  # rho and ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial_cont,
    ald = ald_trial_cont,
    rho = rho_new,
    trt_var = "trt",
    N = 10000,
    marginal_distns = marginals_orig$marginal_dists,
    marginal_params = marginals_orig$marginal_params
  )
  
  expect_equal(
    object = mean(result$EM_cont_1),
    expected = marginals_orig$marginal_params$EM_cont_1$mean,
    tolerance = 0.02
  )
  
  expect_equal(
    object = mean(result$EM_cont_2),
    expected = marginals_orig$marginal_params$EM_cont_2$mean,
    tolerance = 0.02
  )
  
  expect_equal(
    object = mean(result$PF_cont_1),
    expected = marginals_orig$marginal_params$PF_cont_1$mean,
    tolerance = 0.02
  )
  
  expect_equal(
    object = mean(result$PF_cont_2),
    expected = marginals_orig$marginal_params$PF_cont_2$mean,
    tolerance = 0.02
  )
  
  expect_equal(
    object = cor(result$EM_cont_1, result$EM_cont_2),
    expected = rho_new,
    tolerance = 0.1
  )
  expect_equal(
    object = cor(result$PF_cont_1, result$PF_cont_2),
    expected = rho_new,
    tolerance = 0.1
  )
  
  # no ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ald = ald_trial_cont,
    rho = rho_new,
    trt_var = "trt",
    N = 10000,
    marginal_distns = marginals_orig$marginal_dists,
    marginal_params = marginals_orig$marginal_params
  )
  
  expect_equal(
    object = mean(result$EM_cont_1),
    expected = marginals_orig$marginal_params$EM_cont_1$mean,
    tolerance = 0.1
  )
  
  expect_equal(
    object = mean(result$EM_cont_2),
    expected = marginals_orig$marginal_params$EM_cont_2$mean,
    tolerance = 0.1
  )
  
  expect_equal(
    object = mean(result$PF_cont_1),
    expected = marginals_orig$marginal_params$PF_cont_1$mean,
    tolerance = 0.1
  )
  
  expect_equal(
    object = mean(result$PF_cont_2),
    expected = marginals_orig$marginal_params$PF_cont_2$mean,
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$EM_cont_1, result$EM_cont_2),
    expected = rho_new,
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$PF_cont_1, result$PF_cont_2),
    expected = rho_new,
    tolerance = 0.1
  )
})

test_that("simulate_ALD_pseudo_pop continuous via outstandR", {
  ##TODO:
  # rho = rho_new,
  # trt_var = "trt",
  # N = 10000,
  
  form <- y ~ PF_cont_1 + PF_cont_2 + trt*EM_cont_1 + trt*EM_cont_2
  corX <- 0.2
  
  meanX_BC <- c(0.6, 0.6)  # mean of each normally-distributed covariate in BC
  meanX_EM_BC <- c(0.6, 0.6)  # mean of each normally-distributed EM covariate in BC
  sdX_EM <- c(0.4, 0.4)  # standard deviation of each EM covariate
  sdX <- c(0.4, 0.4)  # standard deviation of each covariate (same for AC and BC)
  
  marginals_orig <- list(
    marginal_dists = c(EM_cont_1 = "norm",
                       EM_cont_2 = "norm",
                       PF_cont_1 = "norm",
                       PF_cont_2 = "norm"),
    marginal_params = list(
      EM_cont_1 = list(mean = meanX_EM_BC[1], sd = sdX_EM[1]),
      EM_cont_2 = list(mean = meanX_EM_BC[2], sd = sdX_EM[2]),
      PF_cont_1 = list(mean = meanX_BC[1], sd = sdX[1]),
      PF_cont_2 = list(mean = meanX_BC[2], sd = sdX[2])
    )
  )
  
  res <- outstandR(
    ipd_trial = ipd_trial_cont,
    ald_trial = ald_trial_cont,
    strategy = strategy_gcomp_ml(
      formula = form,
      marginal_distns = marginals_orig$marginal_dists,
      marginal_params = marginals_orig$marginal_params))
  
  # expect_equivalent(res)
})


test_that("simulate directly with binary data without marginals", {
  
  form <- y ~ PF_bin_1 + PF_bin_2 + trt:EM_bin_1 + trt:EM_bin_2
  corX <- 0.2
  
  # providing competing arguments
  
  # no rho ie taken from ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial_bin,
    ald = ald_trial_bin,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial_bin |>
    dplyr::filter(variable == "EM_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial_bin |>
    dplyr::filter(variable == "EM_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial_bin |>
    dplyr::filter(variable == "PF_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial_bin |>
    dplyr::filter(variable == "PF_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$EM_bin_1),
    expected = EM1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$EM_bin_2),
    expected = EM2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_bin_1),
    expected = PF1_mean,
    tolerance = 0.01
  )
  expect_equal(
    object = mean(result$PF_bin_2),
    expected = PF2_mean,
    tolerance = 0.01
  )
  
  # ##TODO:
  # expect_equal(
  #   object = cor(result$EM_bin_1, result$EM_bin_2),
  #   expected = cor(ipd_trial$EM_bin_1, ipd_trial$EM_bin_2),
  #   tolerance = 0.1
  # )
  # 
  # expect_equal(
  #   object = cor(result$PF_bin_1, result$PF_bin_2),
  #   expected = cor(ipd_trial$PF_bin_1, ipd_trial$PF_bin_2),
  #   tolerance = 0.1
  # )
  
  # rho (same) and ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial_bin,
    ald = ald_trial_bin,
    rho = corX,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "EM_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "EM_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "PF_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "PF_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$EM_bin_1),
    expected = EM1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$EM_bin_2),
    expected = EM2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_bin_1),
    expected = PF1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_bin_2),
    expected = PF2_mean,
    tolerance = 0.01
  )
  
  # expect_equal(
  #   object = cor(result$PF_bin_1, result$PF_bin_2),
  #   expected = corX,
  #   tolerance = 0.1
  # )
  # 
  # expect_equal(
  #   object = cor(result$EM_bin_1, result$EM_bin_2),
  #   expected = corX,
  #   tolerance = 0.1
  # )
  
  # rho (different) and ipd
  
  rho_new <- 0.6
  
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial_bin,
    ald = ald_trial_bin,
    rho = rho_new,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "EM_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "EM_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "PF_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "PF_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$EM_bin_1),
    expected = EM1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$EM_bin_2),
    expected = EM2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_bin_1),
    expected = PF1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_bin_2),
    expected = PF2_mean,
    tolerance = 0.01
  )
  
  # expect_equal(
  #   object = cor(result$EM_bin_1, result$EM_bin_2),
  #   expected = rho_new,
  #   tolerance = 0.1
  # )
  # 
  # expect_equal(
  #   object = cor(result$PF_bin_1, result$PF_bin_2),
  #   expected = rho_new,
  #   tolerance = 0.1
  # )
  
  # no ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ald = ald_trial_bin,
    rho = rho_new,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "EM_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "EM_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "PF_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial_bin |> 
    dplyr::filter(variable == "PF_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$EM_bin_1),
    expected = EM1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$EM_bin_2),
    expected = EM2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_bin_1),
    expected = PF1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_bin_2),
    expected = PF2_mean,
    tolerance = 0.01
  )
  
  # expect_equal(
  #   object = cor(result$EM_bin_1, result$EM_bin_2),
  #   expected = rho_new,
  #   tolerance = 0.1
  # )
  # 
  # expect_equal(
  #   object = cor(result$PF_bin_1, result$PF_bin_2),
  #   expected = rho_new,
  #   tolerance = 0.1
  # )
  
})

test_that("simulate directly with binary data with marginals", {
  
  form <- y ~ PF_bin_1 + PF_bin_2 + trt:EM_bin_1 + trt:EM_bin_2
  corX <- 0.2
  
  propX_EM_BC <- c(0.6, 0.6)
  propX_BC <- c(0.6, 0.6)
  
  marginals_orig <- list(
    marginal_dists = c(EM_bin_1 = "binom",
                       EM_bin_2 = "binom",
                       PF_bin_1 = "binom",
                       PF_bin_2 = "binom"),
    marginal_params = list(
      EM_bin_1 = list(prob = propX_EM_BC[1]),
      EM_bin_2 = list(prob = propX_EM_BC[2]),
      PF_bin_1 = list(prob = propX_BC[1]),
      PF_bin_2 = list(prob = propX_BC[2])
    )
  )
  
  rho_new <- 0.6
  
  # no rho ie taken from ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial_bin,
    ald = ald_trial_bin,
    trt_var = "trt",
    N = 100000,
    marginal_distns = marginals_orig$marginal_dists,
    marginal_params = marginals_orig$marginal_params
  )
  
  expect_equal(
    object = mean(result$EM_bin_1),
    expected = marginals_orig$marginal_params$EM_bin_1$prob,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$EM_bin_2),
    expected = marginals_orig$marginal_params$EM_bin_2$prob,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_bin_1),
    expected = marginals_orig$marginal_params$PF_bin_1$prob,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$PF_bin_2),
    expected = marginals_orig$marginal_params$PF_bin_2$prob,
    tolerance = 0.01
  )
  
  ##TODO: why are the binomial correlations wrong?
  
  # expect_equal(
  #   object = cor(result$EM_cont_1, result$EM_bin_2),
  #   expected = cor(ipd_trial$EM_cont_1, ipd_trial$EM_bin_2),
  #   tolerance = 0.1
  # )
  # 
  # expect_equal(
  #   object = cor(result$PF_cont_1, result$PF_bin_2),
  #   expected = cor(ipd_trial$PF_cont_1, ipd_trial$PF_bin_2),
  #   tolerance = 0.1
  # )
  
  # rho and ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial_bin,
    ald = ald_trial_bin,
    rho = rho_new,
    trt_var = "trt",
    N = 10000,
    marginal_distns = marginals_orig$marginal_dists,
    marginal_params = marginals_orig$marginal_params
  )
  
  expect_equal(
    object = mean(result$EM_bin_1),
    expected = marginals_orig$marginal_params$EM_bin_1$prob,
    tolerance = 0.1
  )
  
  expect_equal(
    object = mean(result$EM_bin_2),
    expected = marginals_orig$marginal_params$EM_bin_2$prob,
    tolerance = 0.1
  )
  
  expect_equal(
    object = mean(result$PF_bin_1),
    expected = marginals_orig$marginal_params$PF_bin_1$prob,
    tolerance = 0.1
  )
  
  expect_equal(
    object = mean(result$PF_bin_2),
    expected = marginals_orig$marginal_params$PF_bin_2$prob,
    tolerance = 0.1
  )
  
  # expect_equal(
  #   object = cor(result$EM_bin_1, result$EM_bin_2),
  #   expected = rho_new,
  #   tolerance = 0.1
  # )
  # expect_equal(
  #   object = cor(result$PF_bin_1, result$PF_bin_2),
  #   expected = rho_new,
  #   tolerance = 0.1
  # )
  
  # no ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ald = ald_trial_bin,
    rho = rho_new,
    trt_var = "trt",
    N = 10000,
    marginal_distns = marginals_orig$marginal_dists,
    marginal_params = marginals_orig$marginal_params
  )
  
  expect_equal(
    object = mean(result$EM_bin_1),
    expected = marginals_orig$marginal_params$EM_bin_1$prob,
    tolerance = 0.2
  )
  
  expect_equal(
    object = mean(result$EM_bin_2),
    expected = marginals_orig$marginal_params$EM_bin_2$prob,
    tolerance = 0.02
  )
  
  expect_equal(
    object = mean(result$PF_bin_1),
    expected = marginals_orig$marginal_params$PF_bin_1$prob,
    tolerance = 0.02
  )
  
  expect_equal(
    object = mean(result$PF_bin_2),
    expected = marginals_orig$marginal_params$PF_bin_2$prob,
    tolerance = 0.02
  )
  
  # expect_equal(
  #   object = cor(result$EM_bin_1, result$EM_bin_2),
  #   expected = rho_new,
  #   tolerance = 0.1
  # )
  # 
  # expect_equal(
  #   object = cor(result$PF_bin_1, result$PF_bin_2),
  #   expected = rho_new,
  #   tolerance = 0.1
  # )
  
  # via outstandR()
  
  expect_error(
    outstandR(
      ipd_trial = ipd_trial_bin,
      ald_trial = ald_trial_bin,
      strategy = strategy_maic(formula = form),  # doesnt use pseudo data
      marginal_distns = marginals_orig$marginal_dists,
      marginal_params = marginals_orig$marginal_params),
    regexp = "unused arguments")
  
  # expect_equivalent(res)
})

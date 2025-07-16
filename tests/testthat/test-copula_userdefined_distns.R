#


form <- y ~ PF_cont_1 + PF_cont_2 + trt*EM_cont_1 + trt*EM_cont_2

###########################
# generate continuous data

N <- 200
allocation <- 2 / 3      # active treatment vs. placebo allocation ratio (2:1)
b_trt <- log(0.17)     # conditional effect of active treatment vs. common comparator
b_X <- -log(0.5)       # conditional effect of each prognostic variable
b_EM <- -log(0.67)     # conditional interaction effect of each effect modifier
meanX_AC <- c(0.45, 0.45)       # mean of normally-distributed covariate in AC trial
meanX_BC <- c(0.6, 0.6)         # mean of each normally-distributed covariate in BC
meanX_EM_AC <- c(0.45, 0.45)    # mean of normally-distributed EM covariate in AC trial
meanX_EM_BC <- c(0.6, 0.6)      # mean of each normally-distributed EM covariate in BC
sdX <- c(0.4, 0.4)     # standard deviation of each covariate (same for AC and BC)
sdX_EM <- c(0.4, 0.4)  # standard deviation of each EM covariate
corX <- 0.2            # covariate correlation coefficient
b_0 <- -0.6            # baseline intercept coefficient; fixed value

ipd_trial <- simcovariates::gen_data(
  N = N,
  b_trt = b_trt,
  b_X = b_X,
  b_EM = b_EM,
  b_0 = b_0,
  meanX = meanX_AC,
  sdX = sdX,
  meanX_EM = meanX_EM_AC,
  sdX_EM = sdX_EM,
  corX = corX,
  allocation = allocation,
  family = binomial("logit")
)

ipd_trial$trt <- factor(ipd_trial$trt, labels = c("C", "A"))

BC.IPD <- simcovariates::gen_data(
  N = N,
  b_trt = b_trt,
  b_X = b_X,
  b_EM = b_EM,
  b_0 = b_0,
  meanX = meanX_BC,
  sdX = sdX,
  meanX_EM = meanX_EM_BC,
  sdX_EM = sdX_EM,
  corX = corX,
  allocation = allocation,
  family = binomial("logit")
)

BC.IPD$trt <- factor(BC.IPD$trt, labels = c("C", "B"))

# covariate summary statistics
# assume same between treatments
cov.X <-
  BC.IPD |>
  as.data.frame() |>
  dplyr::select(matches("^(PF|EM)"), trt) |> 
  tidyr::pivot_longer(
    cols = starts_with("PF") | starts_with("EM"),
    names_to = "variable",
    values_to = "value") |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean = mean(value), sd = sd(value)) |>
  tidyr::pivot_longer(
    cols = c("mean", "sd"),
    names_to = "statistic",
    values_to = "value"
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(trt = NA)

# outcome
summary.y <-
  BC.IPD |>
  as.data.frame() |>
  dplyr::select(y, trt) |>
  tidyr::pivot_longer(
    cols = "y",
    names_to = "variable",
    values_to = "value") |>
  dplyr::group_by(variable, trt) |>
  dplyr::summarise(mean = mean(value),
                   sd = sd(value),
                   sum = sum(value)) |>
  tidyr::pivot_longer(
    cols = c("mean", "sd", "sum"),
    names_to = "statistic",
    values_to = "value"
  ) |>
  dplyr::ungroup()

# sample sizes
summary.N <-
  BC.IPD |>
  dplyr::group_by(trt) |>
  dplyr::count(name = "N") |>
  tidyr::pivot_longer(
    cols = "N",
    names_to = "statistic",
    values_to = "value") |>
  dplyr::mutate(variable = NA_character_) |>
  dplyr::select(variable, statistic, value, trt)

ald_trial <- rbind.data.frame(cov.X, summary.y, summary.N)


test_that("simulate_ALD_pseudo_pop directly", {
  ####################
  # without marginals
  
  # providing competing arguments
  
  # no rho ie taken from ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial,
    ald = ald_trial,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial |> 
    dplyr::filter(variable == "EM_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial |> 
    dplyr::filter(variable == "EM_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial |> 
    dplyr::filter(variable == "PF_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial |> 
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
    expected = cor(ipd_trial$EM_cont_1, ipd_trial$EM_cont_2),
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$PF_cont_1, result$PF_cont_2),
    expected = cor(ipd_trial$PF_cont_1, ipd_trial$PF_cont_2),
    tolerance = 0.1
  )
  
  # rho (same) and ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial,
    ald = ald_trial,
    rho = corX,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial |>
    dplyr::filter(variable == "EM_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial |> 
    dplyr::filter(variable == "EM_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial |> 
    dplyr::filter(variable == "PF_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial |> 
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
    ipd = ipd_trial,
    ald = ald_trial,
    rho = rho_new,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial |> dplyr::filter(variable == "EM_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial |> dplyr::filter(variable == "EM_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial |> dplyr::filter(variable == "PF_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial |> dplyr::filter(variable == "PF_cont_2", statistic == "mean") |>
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
    ald = ald_trial,
    rho = rho_new,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial |> 
    dplyr::filter(variable == "EM_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial |>
    dplyr::filter(variable == "EM_cont_2", statistic == "mean") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial |> 
    dplyr::filter(variable == "PF_cont_1", statistic == "mean") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial |> 
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
    ipd = ipd_trial,
    ald = ald_trial,
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
    expected = cor(ipd_trial$EM_cont_1, ipd_trial$EM_cont_2),
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$PF_cont_1, result$PF_cont_2),
    expected = cor(ipd_trial$PF_cont_1, ipd_trial$PF_cont_2),
    tolerance = 0.1
  )
  
  # rho and ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial,
    ald = ald_trial,
    rho = rho_new,
    trt_var = "trt",
    N = 10000,
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
    ald = ald_trial,
    rho = rho_new,
    trt_var = "trt",
    N = 10000,
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
    ipd_trial = ipd_trial,
    ald_trial = ald_trial,
    strategy = strategy_gcomp_ml(
      formula = form,
      marginal_distns = marginals_orig$marginal_dists,
      marginal_params = marginals_orig$marginal_params))
  
  # expect_equivalent(res)
})



#######################
# generate binary data

form <- y ~ PF_bin_1 + PF_bin_2 + trt:EM_bin_1 + trt:EM_bin_2

N <- 200
allocation <- 2 / 3    # active treatment vs. placebo allocation ratio (2:1)
b_trt <- log(0.17)     # conditional effect of active treatment vs. common comparator
b_X <- -log(0.5)       # conditional effect of each prognostic variable
b_EM <- -log(0.67)     # conditional interaction effect of each effect modifier
propX_AC <- c(0.45, 0.45)
propX_BC <- c(0.6, 0.6)  
propX_EM_AC <- c(0.45, 0.45)
propX_EM_BC <- c(0.6, 0.6)  
corX <- 0.2            # covariate correlation coefficient
b_0 <- -0.6            # baseline intercept coefficient; fixed value

ipd_trial <- simcovariates::gen_data(
  N = N,
  b_trt = b_trt,
  b_X_bin = b_X,
  b_EM_bin = b_EM,
  prob_X_bin = propX_AC,
  prob_EM_bin = propX_EM_AC,
  b_0 = b_0,
  corX = corX,
  allocation = allocation,
  family = binomial("logit")
)

ipd_trial$trt <- factor(ipd_trial$trt, labels = c("C", "A"))

BC.IPD <- simcovariates::gen_data(
  N = N,
  b_trt = b_trt,
  b_X_bin = b_X,
  b_EM_bin = b_EM,
  prob_X_bin = propX_BC,
  prob_EM_bin = propX_EM_BC,
  b_0 = b_0,
  corX = corX,
  allocation = allocation,
  family = binomial("logit")
)

BC.IPD$trt <- factor(BC.IPD$trt, labels = c("C", "B"))

# covariate summary statistics
# assume same between treatments
cov.X <-
  BC.IPD |>
  as.data.frame() |>
  dplyr::select(matches("^(PF|EM)"), trt) |> 
  tidyr::pivot_longer(
    cols = starts_with("PF") | starts_with("EM"),
    names_to = "variable",
    values_to = "value") |>
  dplyr::group_by(variable) |>
  dplyr::summarise(prop = mean(value)) |>
  tidyr::pivot_longer(
    cols = "prop",
    names_to = "statistic",
    values_to = "value"
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(trt = NA)

# outcome
summary.y <-
  BC.IPD |>
  as.data.frame() |>
  dplyr::select(y, trt) |>
  tidyr::pivot_longer(
    cols = "y",
    names_to = "variable",
    values_to = "value") |>
  dplyr::group_by(variable, trt) |>
  dplyr::summarise(mean = mean(value),
                   sd = sd(value),
                   sum = sum(value)) |>
  tidyr::pivot_longer(
    cols = c("mean", "sd", "sum"),
    names_to = "statistic",
    values_to = "value"
  ) |>
  dplyr::ungroup()

# sample sizes
summary.N <-
  BC.IPD |>
  dplyr::group_by(trt) |>
  dplyr::count(name = "N") |>
  tidyr::pivot_longer(
    cols = "N",
    names_to = "statistic",
    values_to = "value") |>
  dplyr::mutate(variable = NA_character_) |>
  dplyr::select(variable, statistic, value, trt)

ald_trial <- rbind.data.frame(cov.X, summary.y, summary.N)


test_that("simulate directly with binary data without marginals", {
  
  # providing competing arguments
  
  # no rho ie taken from ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial,
    ald = ald_trial,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial |>
    dplyr::filter(variable == "EM_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial |>
    dplyr::filter(variable == "EM_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial |>
    dplyr::filter(variable == "PF_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial |>
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
  
  ##TODO:
  expect_equal(
    object = cor(result$EM_bin_1, result$EM_bin_2),
    expected = cor(ipd_trial$EM_bin_1, ipd_trial$EM_bin_2),
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$PF_bin_1, result$PF_bin_2),
    expected = cor(ipd_trial$PF_bin_1, ipd_trial$PF_bin_2),
    tolerance = 0.1
  )
  
  # rho (same) and ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial,
    ald = ald_trial,
    rho = corX,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial |> 
    dplyr::filter(variable == "EM_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial |> 
    dplyr::filter(variable == "EM_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial |> 
    dplyr::filter(variable == "PF_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial |> 
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
    ipd = ipd_trial,
    ald = ald_trial,
    rho = rho_new,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial |> 
    dplyr::filter(variable == "EM_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial |> 
    dplyr::filter(variable == "EM_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial |> 
    dplyr::filter(variable == "PF_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial |> 
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
    ald = ald_trial,
    rho = rho_new,
    trt_var = "trt",
    N = 100000)
  
  EM1_mean <- ald_trial |> 
    dplyr::filter(variable == "EM_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  EM2_mean <- ald_trial |> 
    dplyr::filter(variable == "EM_bin_2", statistic == "prop") |>
    dplyr::pull(value)
  
  PF1_mean <- ald_trial |> 
    dplyr::filter(variable == "PF_bin_1", statistic == "prop") |>
    dplyr::pull(value)
  
  PF2_mean <- ald_trial |> 
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
    ipd = ipd_trial,
    ald = ald_trial,
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
    ipd = ipd_trial,
    ald = ald_trial,
    rho = rho_new,
    trt_var = "trt",
    N = 10000,
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
    ald = ald_trial,
    rho = rho_new,
    trt_var = "trt",
    N = 10000,
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
      ipd_trial = ipd_trial,
      ald_trial = ald_trial,
      strategy = strategy_maic(formula = form),  # doesnt use pseudo data
      marginal_distns = marginals_orig$marginal_dists,
      marginal_params = marginals_orig$marginal_params),
    regexp = "unused arguments")
  
  # expect_equivalent(res)
})

#


form <- y ~ X3 + X4 + trt:X1 + trt:X2

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
  tidyr::pivot_longer(cols = "y",
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
  tidyr::pivot_longer(cols = "N",
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
    N = 100000
  )
  
  X1_mean <- ald_trial |> dplyr::filter(variable == "X1", statistic == "mean") |>
    dplyr::pull(value)
  
  X2_mean <- ald_trial |> dplyr::filter(variable == "X2", statistic == "mean") |>
    dplyr::pull(value)
  
  X3_mean <- ald_trial |> dplyr::filter(variable == "X3", statistic == "mean") |>
    dplyr::pull(value)

  X4_mean <- ald_trial |> dplyr::filter(variable == "X4", statistic == "mean") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$X1),
    expected = X1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X2),
    expected = X2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X3),
    expected = X3_mean,
    tolerance = 0.01
  )
  expect_equal(
    object = mean(result$X4),
    expected = X4_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$X3, result$X4),
    expected = cor(ipd_trial$X3, ipd_trial$X4),
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$X1, result$X2),
    expected = cor(ipd_trial$X1, ipd_trial$X2),
    tolerance = 0.1
  )
  
  # rho (same) and ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ipd = ipd_trial,
    ald = ald_trial,
    rho = corX,
    trt_var = "trt",
    N = 100000
  )
  
  X1_mean <- ald_trial |> dplyr::filter(variable == "X1", statistic == "mean") |>
    dplyr::pull(value)
  
  X2_mean <- ald_trial |> dplyr::filter(variable == "X2", statistic == "mean") |>
    dplyr::pull(value)
  
  X3_mean <- ald_trial |> dplyr::filter(variable == "X3", statistic == "mean") |>
    dplyr::pull(value)
  
  X4_mean <- ald_trial |> dplyr::filter(variable == "X4", statistic == "mean") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$X1),
    expected = X1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X2),
    expected = X2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X3),
    expected = X3_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X4),
    expected = X4_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$X3, result$X4),
    expected = corX,
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$X1, result$X2),
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
    N = 100000
  )
  
  X1_mean <- ald_trial |> dplyr::filter(variable == "X1", statistic == "mean") |>
    dplyr::pull(value)
  
  X2_mean <- ald_trial |> dplyr::filter(variable == "X2", statistic == "mean") |>
    dplyr::pull(value)
  
  X3_mean <- ald_trial |> dplyr::filter(variable == "X3", statistic == "mean") |>
    dplyr::pull(value)
  
  X4_mean <- ald_trial |> dplyr::filter(variable == "X4", statistic == "mean") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$X1),
    expected = X1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X2),
    expected = X2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X3),
    expected = X3_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X4),
    expected = X4_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$X3, result$X4),
    expected = rho_new,
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$X1, result$X2),
    expected = rho_new,
    tolerance = 0.1
  )
  
  # no ipd
  result <- simulate_ALD_pseudo_pop(
    formula = form,
    ald = ald_trial,
    rho = rho_new,
    trt_var = "trt",
    N = 100000
  )
  
  X1_mean <- ald_trial |> dplyr::filter(variable == "X1", statistic == "mean") |>
    dplyr::pull(value)
  
  X2_mean <- ald_trial |> dplyr::filter(variable == "X2", statistic == "mean") |>
    dplyr::pull(value)
  
  X3_mean <- ald_trial |> dplyr::filter(variable == "X3", statistic == "mean") |>
    dplyr::pull(value)
  
  X4_mean <- ald_trial |> dplyr::filter(variable == "X4", statistic == "mean") |>
    dplyr::pull(value)
  
  expect_equal(
    object = mean(result$X1),
    expected = X1_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X2),
    expected = X2_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X3),
    expected = X3_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X4),
    expected = X4_mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$X3, result$X4),
    expected = rho_new,
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$X1, result$X2),
    expected = rho_new,
    tolerance = 0.1
  )
  
  #################
  # with marginals
  
  marginals_orig <- list(
    marginal_dists = c(X1 = "norm", X2 = "norm", X3 = "norm", X4 = "norm"),
    marginal_params = list(
      X1 = list(mean = meanX_EM_BC[1], sd = sdX_EM[1]),
      X2 = list(mean = meanX_EM_BC[2], sd = sdX_EM[2]),
      X3 = list(mean = meanX_BC[1], sd = sdX[1]),
      X4 = list(mean = meanX_BC[2], sd = sdX[2])
    )
  )
  
  marginals_new <- list(
    marginal_dists = c(X1 = "norm", X2 = "norm", X3 = "norm", X4 = "norm"),
    marginal_params = list(
      X1 = list(mean = 0, sd = 1),
      X2 = list(mean = 0, sd = 1),
      X3 = list(mean = 0, sd = 1),
      X4 = list(mean = 0, sd = 1)
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
    object = mean(result$X1),
    expected = marginals_orig$marginal_params[[1]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X2),
    expected = marginals_orig$marginal_params[[2]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X3),
    expected = marginals_orig$marginal_params[[3]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X4),
    expected = marginals_orig$marginal_params[[4]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$X3, result$X4),
    expected = cor(ipd_trial$X3, ipd_trial$X4),
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$X1, result$X2),
    expected = cor(ipd_trial$X1, ipd_trial$X2),
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
    object = mean(result$X1),
    expected = marginals_orig$marginal_params[[1]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X2),
    expected = marginals_orig$marginal_params[[2]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X3),
    expected = marginals_orig$marginal_params[[3]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X4),
    expected = marginals_orig$marginal_params[[4]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$X3, result$X4),
    expected = rho_new,
    tolerance = 0.1
  )
  expect_equal(
    object = cor(result$X1, result$X2),
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
    object = mean(result$X1),
    expected = marginals_orig$marginal_params[[1]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X2),
    expected = marginals_orig$marginal_params[[2]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X3),
    expected = marginals_orig$marginal_params[[3]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = mean(result$X4),
    expected = marginals_orig$marginal_params[[4]]$mean,
    tolerance = 0.01
  )
  
  expect_equal(
    object = cor(result$X3, result$X4),
    expected = rho_new,
    tolerance = 0.1
  )
  
  expect_equal(
    object = cor(result$X1, result$X2),
    expected = rho_new,
    tolerance = 0.1
  )
})

#######################
# generate binary data

N <- 200
allocation <- 2 / 3      # active treatment vs. placebo allocation ratio (2:1)
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
  prop_X_bin = propX_AC,
  prop_EM_bin = propX_EM_AC,
  b_0 = b_0,
  corX = corX,
  allocation = allocation,
  family = binomial("logit")
)

BC.IPD <- simcovariates::gen_data(
  N = N,
  b_trt = b_trt,
  b_X_bin = b_X,
  b_EM_bin = b_EM,
  prop_X_bin = propX_BC,
  prop_EM_bin = propX_EM_BC,
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
  tidyr::pivot_longer(cols = "y",
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
  tidyr::pivot_longer(cols = "N",
               names_to = "statistic",
               values_to = "value") |>
  dplyr::mutate(variable = NA_character_) |>
  dplyr::select(variable, statistic, value, trt)

ald_trial <- rbind.data.frame(cov.X, summary.y, summary.N)


test_that("simulate_ALD_pseudo_pop directly with binary data", {
  ##TODO:
})














test_that("simulate_ALD_pseudo_pop via outstandR", {
  ##TODO:
  # outstandR(marginal_distns =  c("norm", "exp"),
  #           marginal_params = list(
  #     list(mean = 55, sd = 8),
  #     list(rate = 0.04)
  #   )
  
  # expect_
})

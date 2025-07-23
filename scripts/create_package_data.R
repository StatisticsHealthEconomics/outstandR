# create package data


library(simcovariates)
library(dplyr)

set.seed(220725)

N <- 200
allocation <- 2/3      # active treatment vs. placebo allocation ratio (2:1)
b_trt <- log(0.17)     # conditional effect of active treatment vs. common comparator
b_0 <- -0.6            # baseline intercept coefficient  ##TODO: fixed value


####################################
# binary outcome, continuous covariates

b_PF <- -log(0.5)       # conditional effect of each prognostic variable
b_EM <- -log(0.67)      # conditional interaction effect of each effect modifier
meanX_AC <- c(0.45, 0.45)       # mean of normally-distributed covariate in AC trial
meanX_BC <- c(0.6, 0.6)         # mean of each normally-distributed covariate in BC
meanX_EM_AC <- c(0.45, 0.45)    # mean of normally-distributed EM covariate in AC trial
meanX_EM_BC <- c(0.6, 0.6)      # mean of each normally-distributed EM covariate in BC
sdX <- c(0.4, 0.4)     # standard deviation of each covariate (same for AC and BC)
sdX_EM <- c(0.4, 0.4)  # standard deviation of each EM covariate
corX <- 0.2            # covariate correlation coefficient  

# y ~ PF1 + PF2 + trt + trt:(EM1 + EM2)

covariate_defns_ipd <- list(
  PF_cont_1 = list(type = continuous(mean = meanX_AC[1], sd = sdX[1]),
                   role = "prognostic"),
  PF_cont_2 = list(type = continuous(mean = meanX_AC[2], sd = sdX[2]),
                   role = "prognostic"),
  EM_cont_1 = list(type = continuous(mean = meanX_EM_AC[1], sd = sdX_EM[1]),
                   role = "effect_modifier"),
  EM_cont_2 = list(type = continuous(mean = meanX_EM_AC[2], sd = sdX_EM[2]),
                   role = "effect_modifier")
)

b_prognostic <- c(PF_cont_1 = b_PF, PF_cont_2 = b_PF)

b_effect_modifier <- c(EM_cont_1 = b_EM, EM_cont_2 = b_EM)

num_normal_covs <- length(covariate_defns_ipd)
cor_matrix <- matrix(corX, num_normal_covs, num_normal_covs)
diag(cor_matrix) <- 1

rownames(cor_matrix) <- c("PF_cont_1", "PF_cont_2", "EM_cont_1", "EM_cont_2")
colnames(cor_matrix) <- c("PF_cont_1", "PF_cont_2", "EM_cont_1", "EM_cont_2")

AC_IPD_binY_contX <- simcovariates::gen_data(
  N = N,
  b_0 = b_0,
  b_trt = b_trt,
  covariate_defns = covariate_defns_ipd,
  b_prognostic = b_prognostic,
  b_effect_modifier = b_effect_modifier,
  cor_matrix = cor_matrix,
  trt_assignment = list(prob_trt1 = allocation),
  family = binomial("logit"))

AC_IPD_binY_contX$trt <- factor(AC_IPD_binY_contX$trt, labels = c("C", "A"))

## aggregate-level data

covariate_defns_ald <- list(
  PF_cont_1 = list(type = continuous(mean = meanX_BC[1], sd = sdX[1]),
                   role = "prognostic"),
  PF_cont_2 = list(type = continuous(mean = meanX_BC[2], sd = sdX[2]),
                   role = "prognostic"),
  EM_cont_1 = list(type = continuous(mean = meanX_EM_BC[1], sd = sdX_EM[1]),
                   role = "effect_modifier"),
  EM_cont_2 = list(type = continuous(mean = meanX_EM_BC[2], sd = sdX_EM[2]),
                   role = "effect_modifier")
)

BC.IPD <- simcovariates::gen_data(
  N = N,
  b_0 = b_0,
  b_trt = b_trt,
  covariate_defns = covariate_defns_ald,
  b_prognostic = b_prognostic,
  b_effect_modifier = b_effect_modifier,
  cor_matrix = cor_matrix,
  trt_assignment = list(prob_trt1 = allocation),
  family = binomial("logit"))

BC.IPD$trt <- factor(BC.IPD$trt, labels = c("C", "B"))

# covariate summary statistics
# assume same between treatments
cov.X <- 
  BC.IPD %>%
  as.data.frame() |> 
  dplyr::select(matches("^(PF|EM)"), trt) |> 
  tidyr::pivot_longer(
    cols = starts_with("PF") | starts_with("EM"),
    names_to = "variable",
    values_to = "value") |>
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    mean = mean(value),
    sd = sd(value)
  ) %>%
  tidyr::pivot_longer(
    cols = c("mean", "sd"),
    names_to = "statistic",
    values_to = "value") %>%
  dplyr::ungroup() |> 
  dplyr::mutate(trt = NA)

# outcome
summary.y <- 
  BC.IPD |> 
  as.data.frame() |> 
  dplyr::select(y, trt) %>%
  tidyr::pivot_longer(cols = "y",
                      names_to = "variable",
                      values_to = "value") %>%
  dplyr::group_by(variable, trt) %>%
  dplyr::summarise(
    mean = mean(value),
    sd = sd(value),
    sum = sum(value)
  ) %>%
  tidyr::pivot_longer(
    cols = c("mean", "sd", "sum"),
    names_to = "statistic",
    values_to = "value") %>%
  dplyr::ungroup()

# sample sizes
summary.N <- 
  BC.IPD |> 
  group_by(trt) |> 
  count(name = "N") |> 
  tidyr::pivot_longer(
    cols = "N",
    names_to = "statistic",
    values_to = "value") |> 
  mutate(variable = NA_character_) |> 
  dplyr::select(variable, statistic, value, trt)

BC_ALD_binY_contX <- rbind.data.frame(cov.X, summary.y, summary.N)

save(AC_IPD_binY_contX, file = "data/AC_IPD_binY_contX.rda")
save(BC_ALD_binY_contX, file = "data/BC_ALD_binY_contX.rda")


#######################################
# continuous outcome, mixed covariates

b_PF <- -log(0.5)       # conditional effect of each prognostic variable
b_EM <- -log(0.67)      # conditional interaction effect of each effect modifier
meanX_AC <- c(0.45, 0.45)       # mean of normally-distributed covariate in AC trial
meanX_BC <- c(0.6, 0.6)         # mean of each normally-distributed covariate in BC
meanX_EM_AC <- c(0.45, 0.45)    # mean of normally-distributed EM covariate in AC trial
meanX_EM_BC <- c(0.6, 0.6)      # mean of each normally-distributed EM covariate in BC
sdX <- c(0.4, 0.4)     # standard deviation of each covariate (same for AC and BC)
sdX_EM <- c(0.4, 0.4)  # standard deviation of each EM covariate
corX <- 0.2            # covariate correlation coefficient  


# y ~ X1 + X3 + X4 + trt + trt:(X2 + X3 + X4)

covariate_defns_ipd <- list(
  X1 = list(type = continuous(mean = meanX_AC[1], sd = sdX[1]),
            role = "prognostic"),
  X2 = list(type = binary(prob = 0.5),
            role = "effect_modifier"),
  X3 = list(type = continuous(mean = meanX_EM_AC[1], sd = sdX_EM[1]),
            role = c("prognostic, effect_modifier")),
  X4 = list(type = binary(prob = 0.2),
            role = c("prognostic", "effect_modifier"))
)

b_prognostic <- c(X1 = b_PF, X3 = b_PF, X4 = b_PF)

b_effect_modifier <- c(X2 = b_EM, X3 = b_EM, X4 = b_EM)

num_normal_covs <- length(covariate_defns_ipd)
cor_matrix <- matrix(corX, num_normal_covs, num_normal_covs)
diag(cor_matrix) <- 1

rownames(cor_matrix) <- c("X1", "X2", "X3", "X4")
colnames(cor_matrix) <- c("X1", "X2", "X3", "X4")

AC_IPD_contY_mixedX <- simcovariates::gen_data(
  N = N,
  b_0 = b_0,
  b_trt = b_trt,
  covariate_defns = covariate_defns_ipd,
  b_prognostic = b_prognostic,
  b_effect_modifier = b_effect_modifier,
  cor_matrix = cor_matrix,
  trt_assignment = list(prob_trt1 = allocation),
  family = gaussian("identity"))

AC_IPD_contY_mixedX$trt <- factor(AC_IPD_contY_mixedX$trt, labels = c("C", "A"))

# aggregate-level data

covariate_defns_ald <- list(
  X1 = list(type = continuous(mean = meanX_BC[1], sd = sdX[1]),
            role = "prognostic"),
  X2 = list(type = binary(prob = 0.7),
            role = "effect_modifier"),
  X3 = list(type = continuous(mean = meanX_EM_BC[1], sd = sdX_EM[1]),
            role = c("prognostic, effect_modifier")),
  X4 = list(type = binary(prob = 0.5),
            role = c("prognostic, effect_modifier"))
)

BC.IPD <- simcovariates::gen_data(
  N = N,
  b_0 = b_0,
  b_trt = b_trt,
  covariate_defns = covariate_defns_ald,
  b_prognostic = b_prognostic,
  b_effect_modifier = b_effect_modifier,
  cor_matrix = cor_matrix,
  trt_assignment = list(prob_trt1 = allocation),
  family = gaussian("identity"))

BC.IPD$trt <- factor(BC.IPD$trt, labels = c("C", "B"))

# covariate summary statistics
# assume same between treatments
variable_types_lookup <- tibble::tibble(
  variable = names(covariate_defns_ald),
  var_type_class = sapply(covariate_defns_ald, function(x) class(x$type)[1])
)

cov.X <- BC.IPD %>%
  as.data.frame() |>
  dplyr::select(matches("^X"), trt) |>
  tidyr::pivot_longer(
    cols = starts_with("X"),
    names_to = "variable",
    values_to = "value"
  ) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    mean = mean(value),
    sd = sd(value),
    .groups = "drop"
  ) %>%
  dplyr::left_join(variable_types_lookup, by = "variable") %>%
  dplyr::mutate(
    is_binary = var_type_class == "bin_var_type"
  ) %>%
  tidyr::pivot_longer(
    cols = c("mean", "sd"),
    names_to = "stat_type",
    values_to = "calculated_value"
  ) %>%
  dplyr::filter(
    (is_binary & stat_type == "mean") | # For binary, only keep the 'mean' (which is the probability)
      (!is_binary & stat_type %in% c("mean", "sd"))
  ) %>%
  dplyr::mutate(
    statistic = ifelse(is_binary & stat_type == "mean", "prop", # Rename mean to proportion for binary
                       ifelse(stat_type == "mean", "mean", "sd"))
  ) %>%
  dplyr::select(variable, statistic, value = calculated_value) %>%
  dplyr::mutate(trt = NA)

# outcome
summary.y <- 
  BC.IPD |> 
  as.data.frame() |> 
  dplyr::select(y, trt) %>%
  tidyr::pivot_longer(cols = "y",
                      names_to = "variable",
                      values_to = "value") %>%
  group_by(variable, trt) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    sum = sum(value)
  ) %>%
  tidyr::pivot_longer(
    cols = c("mean", "sd", "sum"),
    names_to = "statistic",
    values_to = "value") %>%
  ungroup()

# sample sizes
summary.N <- 
  BC.IPD |> 
  group_by(trt) |> 
  count(name = "N") |> 
  tidyr::pivot_longer(
    cols = "N",
    names_to = "statistic",
    values_to = "value") |> 
  mutate(variable = NA_character_) |> 
  dplyr::select(variable, statistic, value, trt)

BC_ALD_contY_mixedX <- rbind.data.frame(cov.X, summary.y, summary.N)

save(AC_IPD_contY_mixedX, file = "data/AC_IPD_contY_mixedX.rda")
save(BC_ALD_contY_mixedX, file = "data/BC_ALD_contY_mixedX.rda")


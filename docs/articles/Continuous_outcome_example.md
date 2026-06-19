# Continuous data example

## Introduction

This is the vignette for performing population adjustment methods with
continuous data, in order to compare marginal treatment effects when
there are cross-trial differences in effect modifiers and limited
patient-level data. We will demonstrate how to apply MAIC, G-computation
with ML, G-computation with Bayesian inference and multiple imputation
marginalisation. The document structure follow the binary data example
vignette which should be referred to for more details.

## Example analysis

First, let us load necessary packages.

``` r

# install.packages("outstandR",
#  repos = c("https://statisticshealtheconomics.r-universe.dev", "https://cloud.r-project.org"))
#
# install.packages("simcovariates",
#  repos = c("https://n8thangreen.r-universe.dev", "https://cloud.r-project.org"))

library(boot)      # non-parametric bootstrap in MAIC and ML G-computation
library(copula)    # simulating BC covariates from Gaussian copula
library(rstanarm)  # fit outcome regression, draw outcomes in Bayesian G-computation
library(tidyr)
library(dplyr)
library(MASS)
library(outstandR)
library(simcovariates)
```

### Data

We first simulate both the IPD and ALD data. See the binary data example
vignette for more details on how this is implemented.

``` r

N <- 200
allocation <- 2/3      # active treatment vs. placebo allocation ratio (2:1)
b_trt <- log(0.17)     # conditional effect of active treatment vs. common comparator
b_PF <- -log(0.5)      # conditional effect of each prognostic variable
b_EM <- -log(0.67)     # conditional interaction effect of each effect modifier
meanX_AC <- c(0.45, 0.45)       # mean of normally-distributed covariate in AC trial
meanX_BC <- c(0.6, 0.6)         # mean of each normally-distributed covariate in BC
meanX_EM_AC <- c(0.45, 0.45)    # mean of normally-distributed EM covariate in AC trial
meanX_EM_BC <- c(0.6, 0.6)      # mean of each normally-distributed EM covariate in BC
sdX <- c(0.4, 0.4)     # standard deviation of each covariate (same for AC and BC)
sdX_EM <- c(0.4, 0.4)  # standard deviation of each EM covariate
corX <- 0.2            # covariate correlation coefficient  
b_0 <- -0.6            # baseline intercept coefficient  ##TODO: fixed value
```

The difference with that example is that we change the `family` argument
in `gen_data()` to `gaussian(link = "identity")`, corresponding to the
continuous data case. The `gen_data()` function is available in the
[simcovariates](https://github.com/n8thangreen/simcovariates) package on
GitHub.

``` r

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

ipd_trial <- simcovariates::gen_data(
  N = N,
  b_0 = b_0,
  b_trt = b_trt,
  covariate_defns = covariate_defns_ipd,
  b_prognostic = b_prognostic,
  b_effect_modifier = b_effect_modifier,
  cor_matrix = cor_matrix,
  trt_assignment = list(prob_trt1 = allocation),
  family = gaussian("identity"))

ipd_trial$trt <- factor(ipd_trial$trt, labels = c("C", "A"))
```

Similarly, for the aggregate data but with the additional summarise step
(see binary data example vignette for code).

This general format of the data sets are in a ‘long’ style consisting of
the following.

#### `ipd_trial`: Individual patient data

- `PF_*`: Patient measurements prognostic factors
- `EM_*`: Patient measurements effect modifiers
- `trt`: Treatment label (factor)
- `y`: Continuous numeric

#### `ald_trial`: Aggregate-level data

- `variable`: Covariate name. In the case of treatment arm sample size
  this is `NA`
- `statistic`: Summary statistic name from mean, standard deviation or
  sum
- `value`: Numerical value of summary statistic
- `trt`: Treatment label. Because we assume a common covariate
  distribution between treatment arms this is `NA`

Our data look like the following.

``` r

head(ipd_trial)
```

There are 4 correlated continuous covariates generated per subject,
simulated from a multivariate normal distribution. Treatment `trt` takes
either new treatment *A* or standard of care / status quo *C*. The ITC
is ‘anchored’ via *C*, the common treatment.

``` r

ald_trial
```

In this case, we have 4 covariate mean and standard deviation values;
and the total, average and sample size for each treatment *B* and *C*.

In the following we will implement for MAIC and G-computation methods to
obtain the *marginal variance* and the *marginal treatment effect*.

## Model fitting in R

The [outstandR](https://StatisticsHealthEconomics.github.io/outstandR/)
package has been written to be easy to use and essential consists of a
single function,
[`outstandR()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/outstandR.md).
This can be used to run all of the different types of model, which we
will call *strategies*. The first two arguments of
[`outstandR()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/outstandR.md)
are the individual and aggregate-level data, respectively.

A `strategy` argument of `outstandR` takes functions called
`strategy_*()`, where the wildcard `*` is replaced by the name of the
particular method required,
e.g. [`strategy_maic()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
for MAIC. Each specific example is provided below.

The formula used in this model, passed as an argument to the strategy
function is

``` r

lin_form <- list(
  outcome_model = as.formula("y ~ PF_cont_1 + PF_cont_2 + trt:EM_cont_1 + trt:EM_cont_2"),
  balance_model = as.formula("~ PF_cont_1 + PF_cont_2 + EM_cont_1 + EM_cont_2")
)
```

### MAIC

As mentioned above, pass the model specific strategy function to the
main
[`outstandR()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/outstandR.md)
function, in this case use
[`strategy_maic()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md).

``` r

outstandR_maic <-
  outstandR(ipd_trial, ald_trial,
            strategy = strategy_maic(
              formula = lin_form,
              family = gaussian(link = "identity")))
```

The returned object is of class `outstandR`.

``` r

outstandR_maic
```

### Bayesian G-computation with MCMC

The difference between Bayesian G-computation and its maximum-likelihood
counterpart is in the estimated distribution of the predicted outcomes.
The Bayesian approach also marginalizes, integrates or standardizes over
the joint posterior distribution of the conditional nuisance parameters
of the outcome regression, as well as the joint covariate distribution.

Pass the
[`strategy_gcomp_bayes()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
strategy function.

``` r

outstandR_gcomp_bayes <-
  outstandR(ipd_trial, ald_trial,
            strategy = strategy_gcomp_bayes(
              formula = lin_form,
              family = gaussian(link = "identity")))
```

``` r

outstandR_gcomp_bayes
```

### Multiple imputation marginalisation

Finally, the strategy function to pass to
[`outstandR()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/outstandR.md)
for multiple imputation marginalisation is
[`strategy_mim()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),

``` r

outstandR_mim <-
  outstandR(ipd_trial, ald_trial,
            strategy = strategy_mim(
              formula = lin_form,
              family = gaussian(link = "identity")))
```

``` r

outstandR_mim
```

### Model comparison

Combine all outputs for relative effects table of all contrasts and
methods.

# COPD Analysis

## Analysis

First, let us load necessary packages.

``` r

library(boot)      # non-parametric bootstrap in MAIC and ML G-computation
library(copula)    # simulating BC covariates from Gaussian copula
library(rstanarm)  # fit outcome regression, draw outcomes in Bayesian G-computation
library(outstandR)
```

### Data

``` r

set.seed(555)

AC.IPD <- read.csv(here::here("raw-data", "AC_IPD.csv"))  # AC patient-level data
BC.ALD <- read.csv(here::here("raw-data", "BC_ALD.csv"))  # BC aggregate-level data
```

This general format of data sets consist of the following.

#### `AC.IPD`: Individual patient data

- `X*`: patient measurements
- `trt`: treatment ID
- `y`: (logical) indicator of whether event was observed

#### `BC.ALD`: Aggregate-level data

- `mean.X*`: mean patient measurement
- `sd.X*`: standard deviation of patient measurement
- `y.*.sum`: total number of events
- `y.*.bar`: total number of events
- `N.*`: total number of individuals

Note that the wildcard `*` here is usually an integer from 1 or the
trial identifier *B*, *C*.

Our data look like the following.

``` r

head(AC.IPD)
#>            X1        X2          X3          X4 trt y
#> 1  0.43734111 0.6747901  0.93001035  0.09165363   1 0
#> 2  0.05643081 0.5987971  0.03557646  0.59954129   1 1
#> 3 -0.08048882 0.6843784  0.93147222 -0.11419716   1 0
#> 4 -0.38580926 0.5716644 -0.32252212  0.02551808   1 0
#> 5  1.00755116 0.8220826  0.92735892  0.84414221   1 1
#> 6  0.19443956 0.2031329  0.34990179  0.15633009   1 0
```

There are 4 correlated continuous covariates generated per subject,
simulated from a multivariate normal distribution.

``` r

BC.ALD
#>     mean.X1   mean.X2   mean.X3   mean.X4     sd.X1     sd.X2     sd.X3
#> 1 0.5908996 0.6414179 0.5856529 0.6023671 0.3863145 0.4033615 0.4076097
#>      sd.X4 y.B.sum y.B.bar N.B y.C.sum y.C.bar N.C
#> 1 0.395132     182   0.455 400     149   0.745 200
```

### Output statistics

## Model fitting in R

### MAIC

The formula used in this model is

``` math
y = X_3 + X_4 + \beta_t X_1 + \beta_t X_2
```
which corresponds to the following `R` `formula` object passed as an
argument to the strategy function.

``` r

lin_form <- as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")
```

``` r

outstandR_maic <- outstandR(AC.IPD, BC.ALD, strategy = strategy_maic(formula = lin_form))
```

The returned object is of class `outstandR`.

``` r

outstandR_maic
```

### STC

``` math
g(\mu_n) = \beta_0 + (\boldsymbol{x}_n - \boldsymbol{\theta}) \beta_1 + (\beta_z + (\boldsymbol{x_n^{EM}} - \boldsymbol{\theta^{EM}}) \boldsymbol{\beta_2}) \; \mbox{I}(z_n=1)
```

``` math
y = X_3 + X_4 + \beta_t(X_1 - \bar{X_1}) + \beta_t(X_2 - \bar{X_2})
```
However,
[`outstandR()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/outstandR.md)
knows how to handle this so we can simply pass the same (uncentred)
formula as before.

``` r

outstandR_stc <- outstandR(AC.IPD, BC.ALD, strategy = strategy_stc(formula = lin_form))
outstandR_stc
```

### Parametric G-computation with maximum-likelihood estimation

``` math
g(\mu_n) = \beta_0 + \boldsymbol{x}_n \boldsymbol{\beta_1} + (\beta_z + \boldsymbol{x_n^{EM}} \boldsymbol{\beta_2}) \; \mbox{I}(z_n = 1)
```

``` math
\hat{\mu}_0 = \int_{x^*} g^{-1} (\hat{\beta}_0 + x^* \hat{\beta}_1 ) p(x^*) \; \text{d}x^*
```

``` math
\hat{\Delta}^{(2)}_{10} = g(\hat{\mu}_1) - g(\hat{\mu}_0)
```

``` r

outstandR_gcomp_ml <- outstandR(AC.IPD, BC.ALD, strategy = strategy_gcomp_ml(formula = lin_form))
outstandR_gcomp_ml
```

### Bayesian G-computation with MCMC

``` math
p(y^*_{^z*} \mid \mathcal{D}_{AC}) = \int_{x^*} p(y^* \mid z^*, x^*, \mathcal{D}_{AC}) p(x^* \mid \mathcal{D}_{AC})\; \text{d}x^*
```
``` math
= \int_{x^*} \int_{\beta} p(y^* \mid z^*, x^*, \beta) p(x^* \mid \beta) p(\beta \mid \mathcal{D}_{AC})\; d\beta \; \text{d}x^*
```

``` r

outstandR_gcomp_bayes <- outstandR(AC.IPD, BC.ALD, strategy = strategy_gcomp_bayes(formula = lin_form))
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.7 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 4000 [  0%]  (Warmup)
#> Chain 1: Iteration:  400 / 4000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  800 / 4000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 1200 / 4000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 1600 / 4000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 2000 / 4000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 2001 / 4000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 2400 / 4000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 2800 / 4000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 3200 / 4000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 3600 / 4000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 4000 / 4000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.174 seconds (Warm-up)
#> Chain 1:                0.18 seconds (Sampling)
#> Chain 1:                0.354 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 9e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 4000 [  0%]  (Warmup)
#> Chain 2: Iteration:  400 / 4000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  800 / 4000 [ 20%]  (Warmup)
#> Chain 2: Iteration: 1200 / 4000 [ 30%]  (Warmup)
#> Chain 2: Iteration: 1600 / 4000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 2000 / 4000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 2001 / 4000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 2400 / 4000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 2800 / 4000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 3200 / 4000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 3600 / 4000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 4000 / 4000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.175 seconds (Warm-up)
#> Chain 2:                0.2 seconds (Sampling)
#> Chain 2:                0.375 seconds (Total)
#> Chain 2:
outstandR_gcomp_bayes
```

### Multiple imputation marginalisation

``` r

outstandR_mim <- outstandR(AC.IPD, BC.ALD, strategy = strategy_mim(formula = lin_form))
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 1.6e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.16 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 4000 [  0%]  (Warmup)
#> Chain 1: Iteration:  400 / 4000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  800 / 4000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 1001 / 4000 [ 25%]  (Sampling)
#> Chain 1: Iteration: 1400 / 4000 [ 35%]  (Sampling)
#> Chain 1: Iteration: 1800 / 4000 [ 45%]  (Sampling)
#> Chain 1: Iteration: 2200 / 4000 [ 55%]  (Sampling)
#> Chain 1: Iteration: 2600 / 4000 [ 65%]  (Sampling)
#> Chain 1: Iteration: 3000 / 4000 [ 75%]  (Sampling)
#> Chain 1: Iteration: 3400 / 4000 [ 85%]  (Sampling)
#> Chain 1: Iteration: 3800 / 4000 [ 95%]  (Sampling)
#> Chain 1: Iteration: 4000 / 4000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.093 seconds (Warm-up)
#> Chain 1:                0.3 seconds (Sampling)
#> Chain 1:                0.393 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.1e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.11 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 4000 [  0%]  (Warmup)
#> Chain 2: Iteration:  400 / 4000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  800 / 4000 [ 20%]  (Warmup)
#> Chain 2: Iteration: 1001 / 4000 [ 25%]  (Sampling)
#> Chain 2: Iteration: 1400 / 4000 [ 35%]  (Sampling)
#> Chain 2: Iteration: 1800 / 4000 [ 45%]  (Sampling)
#> Chain 2: Iteration: 2200 / 4000 [ 55%]  (Sampling)
#> Chain 2: Iteration: 2600 / 4000 [ 65%]  (Sampling)
#> Chain 2: Iteration: 3000 / 4000 [ 75%]  (Sampling)
#> Chain 2: Iteration: 3400 / 4000 [ 85%]  (Sampling)
#> Chain 2: Iteration: 3800 / 4000 [ 95%]  (Sampling)
#> Chain 2: Iteration: 4000 / 4000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.085 seconds (Warm-up)
#> Chain 2:                0.288 seconds (Sampling)
#> Chain 2:                0.373 seconds (Total)
#> Chain 2:
outstandR_mim
```

Combine all outputs

``` r

knitr::kable(
  data.frame(
  `MAIC` = unlist(outstandR_maic$contrasts),
  `STC` = unlist(outstandR_stc$contrasts),
  `Gcomp ML` = unlist(outstandR_gcomp_ml$contrasts),
  `Gcomp Bayes` = unlist(outstandR_gcomp_bayes$contrasts),
  `MIM` = unlist(outstandR_mim$contrasts))
)
```

|     |      MAIC |        STC |   Gcomp.ML | Gcomp.Bayes |        MIM |
|:----|----------:|-----------:|-----------:|------------:|-----------:|
| AB  |  0.019758 |  0.0133028 |  0.0140737 |   0.0119835 |  0.0124093 |
| AC  | -0.270242 | -0.2766972 | -0.2759263 |  -0.2780165 | -0.2775907 |
| BC  | -0.290000 | -0.2900000 | -0.2900000 |  -0.2900000 | -0.2900000 |

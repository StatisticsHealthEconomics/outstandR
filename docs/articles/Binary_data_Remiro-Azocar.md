# Binary data example using data from Remiro-Azocar et al. (2020)

## Introduction

This is a simpler version of the Binary Data Example with Simulated Data
example. See this document for more details and exposition.

## General problem

Consider one *AB* trial, for which the company has IPD, and one *AC*
trial, for which only published aggregate data are available. We wish to
estimate a comparison of the effects of treatments *B* and *C* on an
appropriate scale in some target population *P*, denoted by the
parameter $`d_{BC(P)}`$. We make use of bracketed subscripts to denote a
specific population. Within the *AB* population there are parameters
$`\mu_{A( AB)}`$, $`\mu_{B(AB)}`$ and $`\mu_{C(AB)}`$ representing the
expected outcome on each treatment (including parameters for treatments
not studied in the *AB* trial, e.g. treatment *C*). The *AB* trial
provides estimators $`\bar{Y}_{A(AB)}`$ and $`\bar{Y}_{B(AB)}`$ of
$`\mu_{A( AB)}`$, $`\mu_{B(AB)}`$, respectively, which are the summary
outcomes. It is the same situation for the *AC* trial.

For a suitable scale, for example a logit, or risk difference, we form
estimators $`\Delta_{AB(AB)}`$ and $`\Delta_{AC(AC)}`$ of the trial
level (or marginal) relative treatment effects.

``` math

\Delta_{AB(AB)} = g(\bar{Y}_{B{(AB)}}) - g(\bar{Y}_{A{(AB)}})
```

## Example analysis

First, let us load necessary packages.

``` r
library(boot)      # non-parametric bootstrap in MAIC and ML G-computation
library(copula)    # simulating BC covariates from Gaussian copula
library(rstanarm)  # fit outcome regression, draw outcomes in Bayesian G-computation
library(outstandR)
```

### Data

Next, we load the data to use in the analysis. The data comes from a
simulation study in Remiro‐Azócar A, Heath A, Baio G (2020). We consider
binary outcomes using the log-odds ratio as the measure of effect. The
binary outcome may be response to treatment or the occurrence of an
adverse event. For trials *AC* and *BC*, outcome $`y_n`$ for subject
$`n`$ is simulated from a Bernoulli distribution with probabilities of
success generated from logistic regression.

For the *BC* trial, the individual-level covariates and outcomes are
aggregated to obtain summaries. The continuous covariates are summarized
as means and standard deviations, which would be available to the
analyst in the published study in a table of baseline characteristics in
the RCT publication. The binary outcomes are summarized in an overall
event table. Typically, the published study only provides aggregate
information to the analyst.

``` r
set.seed(555)

ipd_trial <- read.csv(here::here("raw-data", "AC_IPD.csv"))  # AC patient-level data
ald_trial <- read.csv(here::here("raw-data", "BC_ALD.csv"))  # BC aggregate-level data
```

This general format of data sets consist of the following.

#### `ipd_trial`: Individual patient data

- `X*`: patient measurements
- `trt`: treatment ID (integer)
- `y`: (logical) indicator of whether event was observed

#### `ald_trial`: Aggregate-level data

- `mean.X*`: mean patient measurement
- `sd.X*`: standard deviation of patient measurement
- `y.*.sum`: total number of events
- `y.*.bar`: proportion of events
- `N.*`: total number of individuals

Note that the wildcard `*` here is usually an integer from 1 or the
trial identifier *B*, *C*.

Let us label the treatment levels

``` r
ipd_trial$trt <- factor(ipd_trial$trt, labels = c("C", "A"))
```

Our data look like the following.

``` r
head(ipd_trial)
#>            X1        X2          X3          X4 trt y
#> 1  0.43734111 0.6747901  0.93001035  0.09165363   A 0
#> 2  0.05643081 0.5987971  0.03557646  0.59954129   A 1
#> 3 -0.08048882 0.6843784  0.93147222 -0.11419716   A 0
#> 4 -0.38580926 0.5716644 -0.32252212  0.02551808   A 0
#> 5  1.00755116 0.8220826  0.92735892  0.84414221   A 1
#> 6  0.19443956 0.2031329  0.34990179  0.15633009   A 0
```

There are 4 correlated continuous covariates generated per subject,
simulated from a multivariate normal distribution.

``` r
ald_trial
#>     mean.X1   mean.X2   mean.X3   mean.X4     sd.X1     sd.X2     sd.X3
#> 1 0.5908996 0.6414179 0.5856529 0.6023671 0.3863145 0.4033615 0.4076097
#>      sd.X4 y.B.sum y.B.bar N.B y.C.sum y.C.bar N.C
#> 1 0.395132     182   0.455 400     149   0.745 200
```

In this case, we have 4 covariate mean and standard deviation values;
and the event total, average and sample size for each treatment *B*, and
*C*.

## Model fitting in R

The [outstandR](https://StatisticsHealthEconomics.github.io/outstandR)
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

### MAIC

Using the individual level data for *AC* firstly we perform
non-parametric bootstrap of the `maic.boot` function with `R = 1000`
replicates. This function fits treatment coefficient for the marginal
effect for *A* vs *C*. The returned value is an object of class `boot`
from the `{boot}` package. We then calculate the bootstrap mean and
variance in the wrapper function `maic_boot_stats`.

The formula used in this model has all covariates as prognostic variable
and is

``` math

y = X_1 + X_2 + X_3 + X_4 + (\beta_t + \beta_x X_1 + \beta_x X_2) t
```

which corresponds to the following `R` `formula` object passed as an
argument to the strategy function.

``` r
lin_form <- as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")
```

``` r
outstandR_maic <- outstandR(ipd_trial, ald_trial,
                            strategy = strategy_maic(formula = lin_form,
                                                     family = binomial(link = "logit")))
```

The returned object is of class `outstandR`.

``` r
outstandR_maic
#> Object of class 'outstandR' 
#> Model: binomial 
#> Scale: log_odds 
#> Common treatment: C 
#> Individual patient data study: AC 
#> Aggregate level data study: BC 
#> Confidence interval level: 0.95 
#> 
#> Contrasts:
#> 
#> # A tibble: 3 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl>      <dbl>      <dbl>
#> 1 AB            0.101    0.173      -0.715      0.917
#> 2 AC           -1.15     0.137      -1.88      -0.427
#> 3 BC           -1.25     0.0364     -1.63      -0.879
#> 
#> Absolute:
#> 
#> # A tibble: 2 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl> <lgl>      <lgl>     
#> 1 A             0.434   0.00255 NA         NA        
#> 2 C             0.704   0.00329 NA         NA
```

### STC

STC is the conventional outcome regression method. It involves fitting a
regression model of outcome on treatment and covariates to the IPD
plugging-in covariate mean values.

``` r
outstandR_stc <- outstandR(ipd_trial, ald_trial,
                           strategy = strategy_stc(formula = lin_form,
                                                   family = binomial(link = "logit")))
outstandR_stc
#> Object of class 'outstandR' 
#> Model: binomial 
#> Scale: log_odds 
#> Common treatment: C 
#> Individual patient data study: AC 
#> Aggregate level data study: BC 
#> Confidence interval level: 0.95 
#> 
#> Contrasts:
#> 
#> # A tibble: 3 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl>      <dbl>      <dbl>
#> 1 AB           -0.176   NA           NA        NA    
#> 2 AC           -1.43    NA           NA        NA    
#> 3 BC           -1.25     0.0364      -1.63     -0.879
#> 
#> Absolute:
#> 
#> # A tibble: 2 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl> <lgl>      <lgl>     
#> 1 A             0.117        NA NA         NA        
#> 2 C             0.356        NA NA         NA
```

For the last two approaches, we perform G-computation firstly with a
frequentist MLE approach and then a Bayesian approach.

### Parametric G-computation with maximum-likelihood estimation

G-computation marginalizes the conditional estimates by separating the
regression modelling from the estimation of the marginal treatment
effect for *A* versus *C*.

``` r
outstandR_gcomp_ml <- outstandR(ipd_trial, ald_trial,
                                strategy = strategy_gcomp_ml(formula = lin_form,
                                                             family = binomial(link = "logit")))
outstandR_gcomp_ml
#> Object of class 'outstandR' 
#> Model: binomial 
#> Scale: log_odds 
#> Common treatment: C 
#> Individual patient data study: AC 
#> Aggregate level data study: BC 
#> Confidence interval level: 0.95 
#> 
#> Contrasts:
#> 
#> # A tibble: 3 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl>      <dbl>      <dbl>
#> 1 AB            0.164    0.147      -0.586      0.915
#> 2 AC           -1.09     0.110      -1.74      -0.438
#> 3 BC           -1.25     0.0364     -1.63      -0.879
#> 
#> Absolute:
#> 
#> # A tibble: 2 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl> <lgl>      <lgl>     
#> 1 A             0.485   0.00235 NA         NA        
#> 2 C             0.734   0.00251 NA         NA
```

### Bayesian G-computation with MCMC

The difference between Bayesian G-computation and its maximum-likelihood
counterpart is in the estimated distribution of the predicted outcomes.
The Bayesian approach also marginalizes, integrates or standardizes over
the joint posterior distribution of the conditional nuisance parameters
of the outcome regression, as well as the joint covariate distribution.

``` r
outstandR_gcomp_bayes <-
  outstandR(ipd_trial, ald_trial,
            strategy = strategy_gcomp_bayes(formula = lin_form,
                                           family = binomial(link = "logit")))
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 7.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.74 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.285 seconds (Warm-up)
#> Chain 1:                0.327 seconds (Sampling)
#> Chain 1:                0.612 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.6e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.16 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.301 seconds (Warm-up)
#> Chain 2:                0.385 seconds (Sampling)
#> Chain 2:                0.686 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 4.4e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.44 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.32 seconds (Warm-up)
#> Chain 3:                0.328 seconds (Sampling)
#> Chain 3:                0.648 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 1.9e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.19 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.298 seconds (Warm-up)
#> Chain 4:                0.262 seconds (Sampling)
#> Chain 4:                0.56 seconds (Total)
#> Chain 4:
outstandR_gcomp_bayes
#> Object of class 'outstandR' 
#> Model: binomial 
#> Scale: log_odds 
#> Common treatment: C 
#> Individual patient data study: AC 
#> Aggregate level data study: BC 
#> Confidence interval level: 0.95 
#> 
#> Contrasts:
#> 
#> # A tibble: 3 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl>      <dbl>      <dbl>
#> 1 AB            0.151    0.141      -0.585      0.887
#> 2 AC           -1.10     0.105      -1.74      -0.467
#> 3 BC           -1.25     0.0364     -1.63      -0.879
#> 
#> Absolute:
#> 
#> # A tibble: 2 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl> <lgl>      <lgl>     
#> 1 A             0.485   0.00218 NA         NA        
#> 2 C             0.736   0.00258 NA         NA
```

### Multiple imputation marginalisation

Fit the model as before.

``` r
outstandR_mim <-
  outstandR(ipd_trial, ald_trial,
            strategy = strategy_mim(formula = lin_form,
                                    family = binomial(link = "logit")))
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.309 seconds (Warm-up)
#> Chain 1:                0.336 seconds (Sampling)
#> Chain 1:                0.645 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.9e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.19 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.309 seconds (Warm-up)
#> Chain 2:                0.335 seconds (Sampling)
#> Chain 2:                0.644 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.2e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.323 seconds (Warm-up)
#> Chain 3:                0.322 seconds (Sampling)
#> Chain 3:                0.645 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 1.8e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.18 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.322 seconds (Warm-up)
#> Chain 4:                0.279 seconds (Sampling)
#> Chain 4:                0.601 seconds (Total)
#> Chain 4:
outstandR_mim
#> Object of class 'outstandR' 
#> Model: binomial 
#> Scale: log_odds 
#> Common treatment: C 
#> Individual patient data study: AC 
#> Aggregate level data study: BC 
#> Confidence interval level: 0.95 
#> 
#> Contrasts:
#> 
#> # A tibble: 3 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl>      <dbl>      <dbl>
#> 1 AB            0.153    0.131      -0.556      0.861
#> 2 AC           -1.10     0.0942     -1.70      -0.499
#> 3 BC           -1.25     0.0364     -1.63      -0.879
#> 
#> Absolute:
#> 
#> # A tibble: 2 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>      <lgl>    <lgl>     <lgl>      <lgl>     
#> 1 A          NA       NA        NA         NA        
#> 2 C          NA       NA        NA         NA
```

### Model comparison

Combine all outputs for log-odds ratio table of all contrasts and
methods.

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

|              |       MAIC |        STC |   Gcomp.ML | Gcomp.Bayes |        MIM |
|:-------------|-----------:|-----------:|-----------:|------------:|-----------:|
| means.AB     |  0.1007707 | -0.1758708 |  0.1641674 |   0.1513437 |  0.1525650 |
| means.AC     | -1.1518383 | -1.4284798 | -1.0884416 |  -1.1012654 | -1.1000440 |
| means.BC     | -1.2526090 | -1.2526090 | -1.2526090 |  -1.2526090 | -1.2526090 |
| variances.AB |  0.1732497 |         NA |  0.1465980 |   0.1410717 |  0.1305757 |
| variances.AC |  0.1368488 |         NA |  0.1101971 |   0.1046708 |  0.0941748 |
| variances.BC |  0.0364009 |  0.0364009 |  0.0364009 |   0.0364009 |  0.0364009 |
| CI.AB1       | -0.7150305 |         NA | -0.5862659 |  -0.5848092 | -0.5556730 |
| CI.AB2       |  0.9165719 |         NA |  0.9146007 |   0.8874966 |  0.8608031 |
| CI.AC1       | -1.8768893 |         NA | -1.7390702 |  -1.7353698 | -1.7015160 |
| CI.AC2       | -0.4267873 |         NA | -0.4378130 |  -0.4671609 | -0.4985720 |
| CI.BC1       | -1.6265510 | -1.6265510 | -1.6265510 |  -1.6265510 | -1.6265510 |
| CI.BC2       | -0.8786671 | -0.8786671 | -0.8786671 |  -0.8786671 | -0.8786671 |

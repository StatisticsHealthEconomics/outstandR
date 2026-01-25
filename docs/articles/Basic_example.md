# Basic Example

## Introduction

Population adjustment methods such as *matching-adjusted indirect
comparison* (MAIC) are increasingly used to compare marginal treatment
effects when there are cross-trial differences in effect modifiers and
limited patient-level data. MAIC is based on propensity score weighting,
which is sensitive to poor covariate overlap and cannot extrapolate
beyond the observed covariate space. Current outcome regression-based
alternatives can extrapolate but target a conditional treatment effect
that is incompatible in the indirect comparison. When adjusting for
covariates, one must integrate or average the conditional estimate over
the relevant population to recover a compatible marginal treatment
effect.

We propose a marginalization method based on *parametric G-computation*
that can be easily applied where the outcome regression is a generalized
linear model or a Cox model. The approach views the covariate adjustment
regression as a nuisance model and separates its estimation from the
evaluation of the marginal treatment effect of interest. The method can
accommodate a Bayesian statistical framework, which naturally integrates
the analysis into a probabilistic framework. A simulation study provides
proof-of-principle and benchmarks the method’s performance against MAIC
and the conventional outcome regression. Parametric G-computation
achieves more precise and more accurate estimates than MAIC,
particularly when covariate overlap is poor, and yields unbiased
marginal treatment effect estimates under no failures of assumptions.
Furthermore, the marginalized regression-adjusted estimates provide
greater precision and accuracy than the conditional estimates produced
by the conventional outcome regression.

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
outcomes. It the same situation for the *AC* trial.

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
the RCT publication. The binary outcomes are summarized as overall event
table. Typically, the published study only provides aggregate
information to the analyst.

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

In this case, we have 4 covariate mean and standard deviation values;
and the event total, average and sample size for each treatment *B*, and
*C*.

### Output statistics

We will implement for MAIC, STC, and G-computation methods to obtain the
*marginal variance*, defined as

``` math

\frac{1}{n_C} + \frac{1}{n_{\bar{C}}} + \frac{1}{n_B} + \frac{1}{n_{\bar{B}}}
```

and the *marginal treatment effect*, defined as

``` math

\log\left( \frac{n_B/(N_B-n_B)}{n_C/(N_B-n_{B})} \right) = \log(n_B n_{\bar{C}}) - \log(n_C n_{\bar{B}})
```
where $`\bar{C}`$ is the compliment of $`C`$ so
e.g. $`n_{\bar{C}} = N_C - n_c`$.

## Model fitting in R

The [outstandR](https://github.com/StatisticsHealthEconomics/outstandR/)
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
non-parametric bootstrap of the `maic.boot` function with `R` = 1000
replicates. This function fits treatment coefficient for the marginal
effect for *A* vs *C*. The returned value is an object of class `boot`
from the `{boot}` package. We then calculate the bootstrap mean and
variance in the wrapper function `maic_boot_stats`.

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
#> $contrasts
#> $contrasts$AB
#> [1] 0.5008915
#> 
#> $contrasts$AC
#> [1] -1.151838
#> 
#> $contrasts$BC
#> [1] -1.65273
#> 
#> 
#> $variances
#> $variances$AB
#>           [,1]
#> [1,] 0.1732497
#> 
#> $variances$AC
#>           [,1]
#> [1,] 0.1368488
#> 
#> $variances$BC
#> [1] 0.03640091
#> 
#> 
#> $CI
#> $CI$AB
#> [1] -0.3149097  1.3166927
#> 
#> $CI$AC
#> [1] -1.8768893 -0.4267873
#> 
#> $CI$BC
#> [1] -2.026672 -1.278788
#> 
#> 
#> attr(,"CI")
#> [1] 0.95
#> attr(,"class")
#> [1] "outstandR" "list"
```

We see that this is a list object with 3 parts, each containing
statistics between each pair of treatments. These are the mean
contrasts, variances and confidence intervals (CI), respectively. The
default CI is for 95% but can be altered in `outstandR` with the `CI`
argument.

### STC

STC is the conventional outcome regression method. It involves fitting a
regression model of outcome on treatment and covariates to the IPD. IPD
effect modifiers are centred at the mean *BC* values.

``` math

g(\mu_n) = \beta_0 + (\boldsymbol{x}_n - \boldsymbol{\theta}) \beta_1 + (\beta_z + (\boldsymbol{x_n^{EM}} - \boldsymbol{\theta^{EM}}) \boldsymbol{\beta_2}) \; \mbox{I}(z_n=1)
```
where $`\beta_0`$ is the intercept, $`\beta_1`$ are the covariate
coefficients, $`\beta_z`$ and $`\beta_2`$ are the effect modifier
coefficients, $`z_n`$ are the indicator variables of effect alternative
treatment. $`g(\cdot)`$ is the link function e.g. $`\log`$.

As already mentioned, running the STC analysis is almost identical to
the previous analysis but we now use the
[`strategy_stc()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
strategy function instead and a formula with centered covariates.

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
#> $contrasts
#> $contrasts$AB
#>     trt 
#> 0.22425 
#> 
#> $contrasts$AC
#>      trt 
#> -1.42848 
#> 
#> $contrasts$BC
#> [1] -1.65273
#> 
#> 
#> $variances
#> $variances$AB
#> [1] 0.1690646
#> 
#> $variances$AC
#> [1] 0.1326637
#> 
#> $variances$BC
#> [1] 0.03640091
#> 
#> 
#> $CI
#> $CI$AB
#> [1] -0.5816375  1.0301375
#> 
#> $CI$AC
#> [1] -2.1423580 -0.7146016
#> 
#> $CI$BC
#> [1] -2.026672 -1.278788
#> 
#> 
#> attr(,"CI")
#> [1] 0.95
#> attr(,"class")
#> [1] "outstandR" "list"
```

For the last two approaches, we perform G-computation firstly with a
frequentist MLE approach and then a Bayesian approach.

### Parametric G-computation with maximum-likelihood estimation

G-computation marginalizes the conditional estimates by separating the
regression modelling from the estimation of the marginal treatment
effect for *A* versus *C*. First, a regression model of the observed
outcome $`y`$ on the covariates $`x`$ and treatment $`z`$ is fitted to
the *AC* IPD:

``` math

g(\mu_n) = \beta_0 + \boldsymbol{x}_n \boldsymbol{\beta_1} + (\beta_z + \boldsymbol{x_n^{EM}} \boldsymbol{\beta_2}) \; \mbox{I}(z_n = 1)
```

In the context of G-computation, this regression model is often called
the “Q-model.” Having fitted the Q-model, the regression coefficients
are treated as nuisance parameters. The parameters are applied to the
simulated covariates $`x*`$ to predict hypothetical outcomes for each
subject under both possible treatments. Namely, a pair of predicted
outcomes, also called potential outcomes, under *A* and under *C*, is
generated for each subject.

By plugging treatment *C* into the regression fit for every simulated
observation, we predict the marginal outcome mean in the hypothetical
scenario in which all units are under treatment *C*:

``` math

\hat{\mu}_0 = \int_{x^*} g^{-1} (\hat{\beta}_0 + x^* \hat{\beta}_1 ) p(x^*) \; \text{d}x^*
```

To estimate the marginal or population-average treatment effect for *A*
versus *C* in the linear predictor scale, one back-transforms to this
scale the average predictions, taken over all subjects on the natural
outcome scale, and calculates the difference between the average linear
predictions:

``` math

\hat{\Delta}^{(2)}_{10} = g(\hat{\mu}_1) - g(\hat{\mu}_0)
```

``` r
outstandR_gcomp_ml <- outstandR(AC.IPD, BC.ALD, strategy = strategy_gcomp_ml(formula = lin_form))
outstandR_gcomp_ml
#> $contrasts
#> $contrasts$AB
#> [1] 0.5645192
#> 
#> $contrasts$AC
#> [1] -1.088211
#> 
#> $contrasts$BC
#> [1] -1.65273
#> 
#> 
#> $variances
#> $variances$AB
#>           [,1]
#> [1,] 0.1425196
#> 
#> $variances$AC
#>           [,1]
#> [1,] 0.1061187
#> 
#> $variances$BC
#> [1] 0.03640091
#> 
#> 
#> $CI
#> $CI$AB
#> [1] -0.1754018  1.3044402
#> 
#> $CI$AC
#> [1] -1.7266857 -0.4497354
#> 
#> $CI$BC
#> [1] -2.026672 -1.278788
#> 
#> 
#> attr(,"CI")
#> [1] 0.95
#> attr(,"class")
#> [1] "outstandR" "list"
```

### Bayesian G-computation with MCMC

The difference between Bayesian G-computation and its maximum-likelihood
counterpart is in the estimated distribution of the predicted outcomes.
The Bayesian approach also marginalizes, integrates or standardizes over
the joint posterior distribution of the conditional nuisance parameters
of the outcome regression, as well as the joint covariate distribution.

Draw a vector of size $`N^*`$ of predicted outcomes $`y^*_z`$ under each
set intervention $`z^* \in \{0, 1\}`$ from its posterior predictive
distribution under the specific treatment. This is defined as
$`p(y^*_{z^*} \mid \mathcal{D}_{AC}) = \int_{\beta} p(y^*_{z^*} \mid \beta) p(\beta \mid \mathcal{D}_{AC}) d\beta`$
where $`p(\beta \mid \mathcal{D}_{AC})`$ is the posterior distribution
of the outcome regression coefficients $`\beta`$, which encode the
predictor-outcome relationships observed in the *AC* trial IPD. This is
given by:

``` math

p(y^*_{^z*} \mid \mathcal{D}_{AC}) = \int_{x^*} p(y^* \mid z^*, x^*, \mathcal{D}_{AC}) p(x^* \mid \mathcal{D}_{AC})\; \text{d}x^*
```
``` math

= \int_{x^*} \int_{\beta} p(y^* \mid z^*, x^*, \beta) p(x^* \mid \beta) p(\beta \mid \mathcal{D}_{AC})\; d\beta \; \text{d}x^*
```

In practice, the integrals above can be approximated numerically, using
full Bayesian estimation via Markov chain Monte Carlo (MCMC) sampling.

The average, variance and interval estimates of the marginal treatment
effect can be derived empirically from draws of the posterior density.

We can draw a vector of size $`N^*`$ of predicted outcomes $`y^*_z`$
under each set intervention $`z^*`$ from its posterior predictive
distribution under the specific treatment.

``` r
outstandR_gcomp_bayes <- outstandR(AC.IPD, BC.ALD, strategy = strategy_gcomp_bayes(formula = lin_form))
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.38 seconds.
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
#> Chain 1:  Elapsed Time: 0.522 seconds (Warm-up)
#> Chain 1:                0.499 seconds (Sampling)
#> Chain 1:                1.021 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.8e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.18 seconds.
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
#> Chain 2:  Elapsed Time: 0.508 seconds (Warm-up)
#> Chain 2:                0.534 seconds (Sampling)
#> Chain 2:                1.042 seconds (Total)
#> Chain 2:
outstandR_gcomp_bayes
#> $contrasts
#> $contrasts$AB
#> [1] 0.5614387
#> 
#> $contrasts$AC
#> [1] -1.091291
#> 
#> $contrasts$BC
#> [1] -1.65273
#> 
#> 
#> $variances
#> $variances$AB
#> [1] 0.137929
#> 
#> $variances$AC
#> [1] 0.1015281
#> 
#> $variances$BC
#> [1] 0.03640091
#> 
#> 
#> $CI
#> $CI$AB
#> [1] -0.1664684  1.2893457
#> 
#> $CI$AC
#> [1] -1.7158038 -0.4667785
#> 
#> $CI$BC
#> [1] -2.026672 -1.278788
#> 
#> 
#> attr(,"CI")
#> [1] 0.95
#> attr(,"class")
#> [1] "outstandR" "list"
```

### Multiple imputation marginalisation

ref

``` math

equation here
```

``` r
outstandR_mim <- outstandR(AC.IPD, BC.ALD, strategy = strategy_mim(formula = lin_form))
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
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
#> Chain 1:  Elapsed Time: 0.254 seconds (Warm-up)
#> Chain 1:                0.751 seconds (Sampling)
#> Chain 1:                1.005 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.8e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.18 seconds.
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
#> Chain 2:  Elapsed Time: 0.262 seconds (Warm-up)
#> Chain 2:                0.768 seconds (Sampling)
#> Chain 2:                1.03 seconds (Total)
#> Chain 2:
outstandR_mim
#> $contrasts
#> $contrasts$AB
#> [1] 0.5811889
#> 
#> $contrasts$AC
#> [1] -1.071541
#> 
#> $contrasts$BC
#> [1] -1.65273
#> 
#> 
#> $variances
#> $variances$AB
#> [1] 0.1330941
#> 
#> $variances$AC
#> [1] 0.09669322
#> 
#> $variances$BC
#> [1] 0.03640091
#> 
#> 
#> $CI
#> $CI$AB
#> [1] -0.1338465  1.2962243
#> 
#> $CI$AC
#> [1] -1.6810021 -0.4620796
#> 
#> $CI$BC
#> [1] -2.026672 -1.278788
#> 
#> 
#> attr(,"CI")
#> [1] 0.95
#> attr(,"class")
#> [1] "outstandR" "list"
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

|     |       MAIC |      STC |   Gcomp.ML | Gcomp.Bayes |        MIM |
|:----|-----------:|---------:|-----------:|------------:|-----------:|
| AB  |  0.5008915 |  0.22425 |  0.5645192 |   0.5614387 |  0.5811889 |
| AC  | -1.1518383 | -1.42848 | -1.0882106 |  -1.0912911 | -1.0715409 |
| BC  | -1.6527298 | -1.65273 | -1.6527298 |  -1.6527298 | -1.6527298 |

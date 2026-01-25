# Continuous data example

\##TODO

## Example analysis

First, let us load necessary packages.

``` r
library(boot)      # non-parametric bootstrap in MAIC and ML G-computation
library(copula)    # simulating BC covariates from Gaussian copula
library(rstanarm)  # fit outcome regression, draw outcomes in Bayesian G-computation
library(outstandR)
library(simcovariates)
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

Use the
[`gen_data()`](https://rdrr.io/pkg/simcovariates/man/gen_data.html)
function available in the
[simcovariates](https://github.com/n8thangreen/simcovariates) package.

``` r
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(MASS)
#> 
#> Attaching package: 'MASS'
#> The following object is masked from 'package:dplyr':
#> 
#>     select

N <- 200
allocation <- 2/3      # active treatment vs. placebo allocation ratio (2:1)
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
b_0 <- -0.6            # baseline intercept coefficient  ##TODO: fixed value

ipd_trial <- gen_data(N, b_trt, b_X, b_EM, b_0,
                      meanX_AC, sdX, 
                      meanX_EM_AC, sdX_EM, 
                      corX, allocation,
                      family = gaussian(link = "identity"))

ipd_trial$trt <- factor(ipd_trial$trt, labels = c("C", "A"))
```

Similarly, for the aggregate data but with the additional summarise
step.

``` r
BC.IPD <- gen_data(N, b_trt, b_X, b_EM, b_0,
                   meanX_BC, sdX, 
                   meanX_EM_BC, sdX_EM, 
                   corX, allocation,
                   family = gaussian(link = "identity"))

cov.X <- BC.IPD %>%
  summarise(across(starts_with("X"),
                   list(mean = mean, sd = sd),
                   .names = "{fn}.{col}"))

out.B <- dplyr::filter(BC.IPD, trt == 1) %>%
  summarise(y.B.sum = sum(y),
            y.B.bar = mean(y),
            y.B.sd = sd(y),
            N.B = n())

out.C <- dplyr::filter(BC.IPD, trt == 0) %>%
  summarise(y.C.sum = sum(y),
            y.C.bar = mean(y),
            y.C.sd = sd(y),
            N.C = n())

ald_trial <- cbind.data.frame(cov.X, out.B, out.C)
```

This general format of data sets consist of the following.

#### `ipd_trial`: Individual patient data

- `X*`: patient measurements
- `trt`: treatment ID (integer)
- `y`: continuous outcome measure

#### `ald_trial`: Aggregate-level data

- `mean.X*`: mean patient measurement
- `sd.X*`: standard deviation of patient measurement
- `mean.y`:
- `sd.y`:

Note that the wildcard `*` here is usually an integer from 1 or the
trial identifier *B*, *C*.

Our data look like the following.

``` r
head(ipd_trial)
#>           X1         X2         X3         X4 trt          y
#> 1 0.42066874  1.0957407 0.37118099  1.3291540   A  0.8949999
#> 2 0.51227329  0.9079984 0.30560144 -0.1842358   A -1.7550940
#> 3 1.00347612  0.8136847 1.25054503  1.1986315   A -0.6521005
#> 4 0.05641685  0.5318140 0.54799796  0.6694090   A -0.6260983
#> 5 0.42997845  0.5274304 0.09140568  0.1222184   A -1.2882834
#> 6 0.06916082 -0.3125090 0.99150666 -0.1102693   A -3.0170689
```

There are 4 correlated continuous covariates generated per subject,
simulated from a multivariate normal distribution.

``` r
ald_trial
#>     mean.X1     sd.X1   mean.X2    sd.X2   mean.X3     sd.X3  mean.X4     sd.X4
#> 1 0.6081015 0.3865186 0.6462269 0.412924 0.5903258 0.3842429 0.635176 0.3831984
#>     y.B.sum   y.B.bar   y.B.sd N.B  y.C.sum   y.C.bar   y.C.sd N.C
#> 1 -135.0634 -1.015514 1.159787 133 12.38076 0.1847875 1.019683  67
```

In this case, we have 4 covariate mean and standard deviation values;
and the event total, average and sample size for each treatment *B*, and
*C*.

### Output statistics

We will implement for MAIC, STC, and G-computation methods to obtain the
*marginal variance*, defined as

``` math
```

and the *marginal treatment effect*, defined as the log-odds ratio,

``` math
```

where $`\bar{C}`$ is the compliment of $`C`$ so e.g.
$`n_{\bar{C}} = N_C - n_c`$.

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
outstandR_maic <-
  outstandR(ipd_trial, ald_trial,
            strategy = strategy_maic(
              formula = lin_form,
              family = gaussian(link = "identity")))
```

The returned object is of class `outstandR`.

``` r
outstandR_maic
#> Object of class 'outstandR' 
#> Model: gaussian 
#> Scale: risk_difference 
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
#> 1 AB           -0.339    0.0533     -0.791      0.113
#> 2 AC           -1.54     0.0276     -1.87      -1.21 
#> 3 BC           -1.20     0.0256     -1.51      -0.887
#> 
#> Absolute:
#> 
#> # A tibble: 2 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl> <lgl>      <lgl>     
#> 1 A            -1.14     0.0112 NA         NA        
#> 2 C             0.396    0.0179 NA         NA
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
outstandR_stc <-
  outstandR(ipd_trial, ald_trial,
            strategy = strategy_stc(
              formula = lin_form,
              family = gaussian(link = "identity")))
outstandR_stc
#> Object of class 'outstandR' 
#> Model: gaussian 
#> Scale: risk_difference 
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
#> 1 AB           -0.363   NA           NA        NA    
#> 2 AC           -1.56    NA           NA        NA    
#> 3 BC           -1.20     0.0256      -1.51     -0.887
#> 
#> Absolute:
#> 
#> # A tibble: 2 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl> <lgl>      <lgl>     
#> 1 A           -1.47          NA NA         NA        
#> 2 C            0.0930        NA NA         NA
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
outstandR_gcomp_ml <-
  outstandR(ipd_trial, ald_trial,
            strategy = strategy_gcomp_ml(formula = lin_form,
                                         family = gaussian(link = "identity")))
outstandR_gcomp_ml
#> Object of class 'outstandR' 
#> Model: gaussian 
#> Scale: risk_difference 
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
#> 1 AB           -0.348    0.0596     -0.826      0.130
#> 2 AC           -1.55     0.0339     -1.91      -1.19 
#> 3 BC           -1.20     0.0256     -1.51      -0.887
#> 
#> Absolute:
#> 
#> # A tibble: 2 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl> <lgl>      <lgl>     
#> 1 A            -1.10     0.0113 NA         NA        
#> 2 C             0.444    0.0241 NA         NA
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
outstandR_gcomp_bayes <-
  outstandR(ipd_trial, ald_trial,
            strategy = strategy_gcomp_bayes(
              formula = lin_form,
              family = gaussian(link = "identity")))
```

``` r
outstandR_gcomp_bayes
#> Object of class 'outstandR' 
#> Model: gaussian 
#> Scale: risk_difference 
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
#> 1 AB           -0.335    0.0560     -0.799      0.129
#> 2 AC           -1.54     0.0304     -1.88      -1.19 
#> 3 BC           -1.20     0.0256     -1.51      -0.887
#> 
#> Absolute:
#> 
#> # A tibble: 2 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl> <lgl>      <lgl>     
#> 1 A            -1.10     0.0113 NA         NA        
#> 2 C             0.433    0.0206 NA         NA
```

### Multiple imputation marginalisation

ref

``` math

equation here
```

``` r
outstandR_mim <-
  outstandR(ipd_trial, ald_trial,
            strategy = strategy_mim(
              formula = lin_form,
              family = gaussian(link = "identity")))
```

``` r
outstandR_mim
#> Object of class 'outstandR' 
#> Model: gaussian 
#> Scale: risk_difference 
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
#> 1 AB           -0.335    0.0538     -0.789      0.120
#> 2 AC           -1.53     0.0281     -1.86      -1.21 
#> 3 BC           -1.20     0.0256     -1.51      -0.887
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

Combine all outputs for relative effects table of all contrasts and
methods.

``` r
knitr::kable(
  data.frame(
    `MAIC` = unlist(outstandR_maic$contrasts$means),
    `STC` = unlist(outstandR_stc$contrasts$means),
    `Gcomp ML` = unlist(outstandR_gcomp_ml$contrasts$means),
    `Gcomp Bayes` = unlist(outstandR_gcomp_bayes$contrasts$means),
    `MIM` = unlist(outstandR_mim$contrasts$means))
  |> 
    round(2))
```

|     |  MAIC |   STC | Gcomp.ML | Gcomp.Bayes |   MIM |
|:----|------:|------:|---------:|------------:|------:|
| AB  | -0.34 | -0.36 |    -0.35 |       -0.34 | -0.33 |
| AC  | -1.54 | -1.56 |    -1.55 |       -1.54 | -1.53 |
| BC  | -1.20 | -1.20 |    -1.20 |       -1.20 | -1.20 |

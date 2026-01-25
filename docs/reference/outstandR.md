# Calculate the difference between treatments using all evidence

This is the main, top-level wrapper for `{outstandR}`. Methods taken
from (Remiro‐Azócar et al. 2022) .

## Usage

``` r
outstandR(
  ipd_trial,
  ald_trial,
  strategy,
  ref_trt = NA,
  CI = 0.95,
  scale = NULL,
  var_method = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- ipd_trial:

  Individual-level patient data. For example, suppose between studies
  *A* and *C*. In a long format and must contain a treatment column and
  outcome column consistent with the formula object. The labels in the
  treatment are used internally so there must be a common treatment with
  the aggregate-level data trial.

- ald_trial:

  Aggregate-level data. For example, suppose between studies *B* and
  *C*. The column names are

  - `variable`: Covariate name. In the case of treatment arm sample size
    this is `NA`,

  - `statistic`: Summary statistic name from "mean", standard deviation
    "sd", probability "prop", or "sum",

  - `value`: Numerical value of summary statistic,

  - `trt`: Treatment label. Because we assume a common covariate
    distribution between treatment arms this is `NA`.

- strategy:

  Computation strategy function. These can be
  [`strategy_maic()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),
  [`strategy_stc()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),
  [`strategy_gcomp_ml()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  and
  [`strategy_gcomp_bayes()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md).

- ref_trt:

  Reference / common / anchoring treatment name.

- CI:

  Confidence interval level; between 0,1 with default 0.95.

- scale:

  Relative treatment effect scale. If `NULL`, the scale is automatically
  determined from the model. Choose from "log-odds",
  "log_relative_risk", "risk_difference", "delta_z", "mean_difference",
  "rate_difference" depending on the data type.

- var_method:

  Variance estimation method.

- seed:

  Random seed.

- ...:

  Additional arguments. Currently, can pass named arguments to
  [`rstanarm::stan_glm()`](https://mc-stan.org/rstanarm/reference/stan_glm.html)
  via
  [`strategy_gcomp_bayes()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md).

## Value

List of length 11 of statistics as a `outstandR` class object.
Containing statistics between each pair of treatments. These are the
mean, variances and confidence intervals, for contrasts and absolute
values.

## References

Remiro‐Azócar A, Heath A, Baio G (2022). “Parametric G‐computation for
compatible indirect treatment comparisons with limited individual
patient data.” *Res. Synth. Methods*, 1–31. ISSN 1759-2879,
[doi:10.1002/jrsm.1565](https://doi.org/10.1002/jrsm.1565) , 2108.12208.

## See also

[`strategy_maic()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
[`strategy_stc()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
[`strategy_gcomp_ml()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
[`strategy_gcomp_bayes()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)

## Examples

``` r
data(AC_IPD_binY_contX)  # A vs C individual patient-level data
data(BC_ALD_binY_contX)  # B vs C aggregate-level data

# linear formula
lin_form <- as.formula("y ~ PF_cont_1 + PF_cont_2 + trt*EM_cont_1 + trt*EM_cont_2")
                                
# sampling values of additional arguments picked for speed
# select appropriate to specific analysis

# matching-adjusted indirect comparison
outstandR_maic <- outstandR(
  AC_IPD_binY_contX, BC_ALD_binY_contX,
  strategy = strategy_maic(formula = lin_form, n_boot = 100))

# simulated treatment comparison
outstandR_stc <- outstandR(
  AC_IPD_binY_contX, BC_ALD_binY_contX,
  strategy = strategy_stc(lin_form))

# \donttest{
# G-computation with maximum likelihood
outstandR_gcomp_ml <- outstandR(
  AC_IPD_binY_contX, BC_ALD_binY_contX,
  strategy = strategy_gcomp_ml(lin_form, n_boot = 100, N =100))

# G-computation with Bayesian inference
outstandR_gcomp_bayes <- outstandR(
  AC_IPD_binY_contX, BC_ALD_binY_contX,
  strategy = strategy_gcomp_bayes(lin_form),
  chains = 1, iter = 1000, warmup = 20)

# Multiple imputation marginalization
outstandR_mim <- outstandR(
  AC_IPD_binY_contX, BC_ALD_binY_contX,
  strategy = strategy_mim(lin_form,
                          N = 100), # size of pseudo-population
  chains = 1, iter = 1000, warmup = 20)
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 1.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.15 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 3
#> Chain 1:            adapt_window = 15
#> Chain 1:            term_buffer = 2
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration:  21 / 1000 [  2%]  (Sampling)
#> Chain 1: Iteration: 120 / 1000 [ 12%]  (Sampling)
#> Chain 1: Iteration: 220 / 1000 [ 22%]  (Sampling)
#> Chain 1: Iteration: 320 / 1000 [ 32%]  (Sampling)
#> Chain 1: Iteration: 420 / 1000 [ 42%]  (Sampling)
#> Chain 1: Iteration: 520 / 1000 [ 52%]  (Sampling)
#> Chain 1: Iteration: 620 / 1000 [ 62%]  (Sampling)
#> Chain 1: Iteration: 720 / 1000 [ 72%]  (Sampling)
#> Chain 1: Iteration: 820 / 1000 [ 82%]  (Sampling)
#> Chain 1: Iteration: 920 / 1000 [ 92%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0 seconds (Warm-up)
#> Chain 1:                0.084 seconds (Sampling)
#> Chain 1:                0.084 seconds (Total)
#> Chain 1: 
#> Warning: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
# }
```

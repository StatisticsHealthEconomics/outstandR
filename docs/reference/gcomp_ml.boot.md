# Bootstrap for G-computation via Maximum Likelihood

For use with a bootstrapping function (e.g.,
[`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html)). On each
bootstrap sample of the data, it calculates a relative treatment effect
(e.g., log odds ratio, log relative risk, or risk difference) using
G-computation with maximum likelihood.

## Usage

``` r
gcomp_ml.boot(
  data,
  indices,
  R,
  outcome_model = NULL,
  family,
  trt_var,
  ref_trt = NA,
  comp_trt = NA,
  rho = NA,
  N = 1000,
  marginal_distns = NA,
  marginal_params = NA,
  ald
)
```

## Arguments

- data:

  A data frame containing the original individual participant data
  (IPD).

- indices:

  A vector of indices supplied by the bootstrapping function, used to
  resample `data`.

- family:

  A 'family' object specifying the distribution and link function (e.g.,
  'binomial'). See stats::family() for more details.

- rho:

  A named square matrix specifying the correlation between covariates
  for synthetic data generation. Defaults to `NA`, assuming
  independence.

- N:

  Synthetic sample size for G-computation

- marginal_distns, marginal_params:

  Marginal distributions and parameters

- ald:

  A data frame of aggregate-level data providing covariate
  distributions.

## Value

A single numeric value representing the relative treatment effect

## See also

[`strategy_gcomp_ml()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)

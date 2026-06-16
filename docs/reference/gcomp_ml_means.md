# G-computation maximum likelihood mean outcomes by arm

G-computation maximum likelihood mean outcomes by arm

## Usage

``` r
gcomp_ml_means(
  outcome_model,
  family,
  ipd,
  ald,
  trt_var,
  rho = NA,
  N = 1000,
  ref_trt,
  comp_trt,
  marginal_distns = NA,
  marginal_params = NA
)
```

## Arguments

- family:

  A 'family' object specifying the distribution and link function (e.g.,
  'binomial'). See stats::family() for more details.

- ipd:

  Individual-level patient data. Dataframe with one row per patient with
  outcome, treatment and covariate columns.

- ald:

  Aggregate-level data. Long format summary statistics for each
  covariate and treatment outcomes. We assume a common distribution for
  each treatment arm.

- rho:

  A named square matrix of covariate correlations; default NA.

- N:

  Synthetic sample size for g-computation

- marginal_distns, marginal_params:

  Marginal distributions and parameters

## Value

A list containing:

- `stats`: Named vector of marginal mean probabilities

- `model`: The fitted glm object

## See also

[`strategy_gcomp_ml()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
[`gcomp_ml.boot()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/gcomp_ml.boot.md)

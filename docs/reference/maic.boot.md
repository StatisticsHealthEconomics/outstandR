# MAIC bootstrap sample

Matching-adjusted indirect comparison bootstrap sampling.

## Usage

``` r
maic.boot(
  ipd,
  indices = 1:nrow(ipd),
  formula,
  family,
  ald,
  trt_var,
  hat_w = NULL
)
```

## Arguments

- ipd:

  Individual-level patient data. Dataframe with one row per patient with
  outcome, treatment and covariate columns.

- indices:

  Vector of indices, same length as original, which define the bootstrap
  sample

- formula:

  Linear regression `formula` object. Prognostic factors (PF) are main
  effects and effect modifiers (EM) are interactions with the treatment
  variable, e.g., y ~ X1 + trt + trt:X2. For covariates as both PF and
  EM use `*` syntax.

- family:

  A 'family' object specifying the distribution and link function (e.g.,
  'binomial'). See stats::family() for more details.

- ald:

  Aggregate-level data. Long format summary statistics for each
  covariate and treatment outcomes. We assume a common distribution for
  each treatment arm.

- hat_w:

  MAIC weights; default `NULL` which calls
  [`maic_weights()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/maic_weights.md)

## Value

Vector of fitted probabilities for treatments *A* and *C*

## See also

[`calc_IPD_stats.maic()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calc_IPD_stats.md)

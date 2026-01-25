# Simulate Aggregate-Level Data Pseudo-Population

Generates a synthetic cohort using a normal copula based on
aggregate-level data.

## Usage

``` r
simulate_ALD_pseudo_pop(
  formula,
  ipd = NULL,
  ald = NULL,
  trt_var,
  rho = NA,
  N = 1000,
  marginal_distns = NA,
  marginal_params = NA,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- formula:

  Linear regression `formula` object. Prognostic factors (PF) are main
  effects and effect modifiers (EM) are interactions with the treatment
  variable, e.g., y ~ X1 + trt + trt:X2. For covariates as both PF and
  EM use `*` syntax.

- ipd:

  Individual-level patient data. Dataframe with one row per patient with
  outcome, treatment and covariate columns.

- ald:

  Aggregate-level data. Long format summary statistics for each
  covariate and treatment outcomes. We assume a common distribution for
  each treatment arm.

- rho:

  A named square matrix of covariate correlations or single value;
  default NA takes from IPD.

- N:

  Sample size for the synthetic cohort. Default is 1000.

- marginal_distns:

  Marginal distributions names; vector default NA. Available
  distributions are given in stats::Distributions. See
  [`copula::Mvdc()`](https://rdrr.io/pkg/copula/man/Mvdc.html) for
  details

- marginal_params:

  Marginal distributions parameters; named list of lists, default NA.
  See [`copula::Mvdc()`](https://rdrr.io/pkg/copula/man/Mvdc.html) for
  details

- seed:

  Random seed

- verbose:

  Default `FALSE`

## Value

A data frame representing the synthetic pseudo-population.

# Multiple imputation marginalization (MIM)

Multiple imputation marginalization (MIM)

## Usage

``` r
calc_mim(strategy, analysis_params, ...)
```

## Arguments

- strategy:

  An object of class `strategy` created by functions such as
  [`strategy_maic()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),
  [`strategy_stc()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),
  or
  [`strategy_mim()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md).
  Contains modelling details like the formula and family.

- ...:

  Additional argument to pass to Stan model

## Value

A list containing:

- `means`: A list containing vectors of posterior means (one per
  synthesis `n_imp`):

  - `A`: Comparator means.

  - `C`: Reference means.

- `model`: A list containing:

  - `fit`: The first-stage
    [`rstanarm::stan_glm()`](https://mc-stan.org/rstanarm/reference/stan_glm.html)
    object.

  - `hats.v`: Vector of variance point estimates for each synthesis.

  - `n_imp`: Number of posterior prediction draws (syntheses).

  - `rho`, `N`, `stan_args`: Strategy and model parameters.

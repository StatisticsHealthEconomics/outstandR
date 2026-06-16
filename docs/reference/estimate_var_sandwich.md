# Estimate Variance Sandwich Estimator

Computes the robust (sandwich) variance estimator for the treatment
effect.

## Usage

``` r
estimate_var_sandwich(strategy, analysis_params, ...)
```

## Arguments

- strategy:

  An object of class `strategy` created by functions such as
  [`strategy_maic()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),
  [`strategy_stc()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),
  or
  [`strategy_mim()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md).
  Contains modelling details like the formula and family.

- analysis_params:

  List of analysis parameters (ipd, ald, etc.)

- ...:

  Additional arguments

## Value

Numeric variance estimate for the treatment contrast

# Input data validator

Input data validator

## Usage

``` r
validate_outstandr(ipd_trial, ald_trial, strategy, CI, scale)
```

## Arguments

- ipd_trial:

  Individual patient data

- ald_trial:

  Aggregate level data

- strategy:

  An object of class `strategy` created by functions such as
  [`strategy_maic()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),
  [`strategy_stc()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),
  or
  [`strategy_mim()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md).
  Contains modelling details like the formula and family.

- CI:

  Confidence interval

- scale:

  Outcome scale

## Value

No return value, called for side effects

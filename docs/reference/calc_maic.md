# Calculate MAIC

Calculate MAIC

## Usage

``` r
calc_maic(strategy, analysis_params)
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

  List of analysis parameters. Must contain `ipd` (individual patient
  data) and `ald` (aggregated lead data).

## Value

A list containing:

- `means`: A list containing:

  - `A`: Bootstrap estimates for comparator treatment group "A".

  - `C`: Bootstrap estimates for reference treatment group "C".

- `model`: A list containing model diagnostics derived from the original
  data:

  - `weights`: Vector of calculated weights for the patients in `ipd`.

  - `ESS`: The Effective Sample Size.

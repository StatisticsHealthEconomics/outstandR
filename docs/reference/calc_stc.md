# Calculate simulated treatment comparison statistics

Calculate simulated treatment comparison statistics

## Usage

``` r
calc_stc(strategy, analysis_params, ...)
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
  data).

- ...:

  Additional arguments.

## Value

A list containing:

- `means`: A list containing:

  - `A`: Mean for comparator treatment group "A".

  - `C`: Mean for reference treatment group "C".

- `model`: The fitted [`stats::glm()`](https://rdrr.io/r/stats/glm.html)
  object.

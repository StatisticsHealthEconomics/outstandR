# Aggregate-level data mean and variance statistics

Computes the mean and variance of marginal treatment effects for
aggregate-level trial data.

## Usage

``` r
calc_ALD_stats(strategy, analysis_params)
```

## Arguments

- strategy:

  A list containing the strategy details, including the family
  distribution.

- analysis_params:

  A list containing:

  - `ald` Aggregate-level trial data

  - `ref_trt` Treatment labels reference (common; e.g. placebo)

  - `ald_comp` Treatment labels comparator

  - `scale` A scaling parameter for the calculation. From "log_odds",
    "risk_difference", "log_relative_risk".

## Value

A list containing:

- `mean`:

  The marginal treatment effect mean.

- `var`:

  The marginal treatment effect variance.

## See also

[`marginal_treatment_effect()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/marginal_treatment_effect.md),
[`marginal_variance()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/marginal_variance.md)

## Examples

``` r
strategy <- list(family = list(family = "binomial"))  # basic version

ald <- data.frame(trt = c("B","C","B","C"),
                  variable = c(NA, NA, "y", "y"),
                  statistic = c("N", "N", "sum", "sum"),
                  value = c(100, 100, 50, 60))

calc_ALD_stats(strategy = strategy,
               list(ald = ald,
                    ref_trt = "C",
                    ald_comp = "B",
                    scale = "log_odds"))
#> $mean
#> [1] -0.4054651
#> 
#> $var
#> [1] 0.08166667
#> 
```

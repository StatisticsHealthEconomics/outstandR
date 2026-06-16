# Aggregate-level data mean and variance statistics

Computes the mean and variance of marginal treatment effects for
aggregate-level trial data.

## Usage

``` r
ALD_stats(strategy, ald, treatments = list("B", "C"), scale)
```

## Arguments

- strategy:

  A list containing the strategy details, including the family
  distribution.

- ald:

  Aggregate-level trial data

- treatments:

  Treatment labels list; default `B`, `C` (common; e.g. placebo)

- scale:

  A scaling parameter for the calculation.

## Value

A list containing:

- mean:

  The marginal treatment effect mean.

- var:

  The marginal treatment effect variance.

## See also

[`marginal_treatment_effect()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/marginal_treatment_effect.md),
[`marginal_variance()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/marginal_variance.md)

## Examples

``` r
if (FALSE) { # \dontrun{
strategy <- list(family = list(family = "binomial"))  # basic version
ald <- data.frame(trial = 1:5, n_B = c(10, 20, 15, 30, 25), n_C = c(12, 18, 20, 25, 22))
ALD_stats(strategy, ald, treatments = list("B", "C"), scale = "log")
} # }
```

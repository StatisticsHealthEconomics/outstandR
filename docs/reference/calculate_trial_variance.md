# Calculate trial variance

Computes the variance of treatment effects for a trial based on the
specified family distribution.

## Usage

``` r
calculate_trial_variance(ald, tid, effect, family)
```

## Arguments

- ald:

  Aggregate-level data. Data frame.

- tid:

  Treatment identifier used to extract relevant columns from `ald`.

- effect:

  A character string specifying the effect scale (e.g., "log_odds",
  "risk_difference").

- family:

  A character string specifying the model family (e.g., "binomial",
  "gaussian").

## Value

Numeric computed variance of treatment effects.

## Examples

``` r
ald <- data.frame(trt = c("B","C","B","C"),
                  variable = c(NA, NA, "y", "y"),
                  statistic = c("N", "N", "sum", "sum"),
                  value = c(100, 100, 50, 60))
                  
calculate_trial_variance(
  ald, tid = "B", effect = "log_odds", family = "binomial")
#> [1] 0.04
  
```

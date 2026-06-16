# Marginal effect variance using the delta method

Computes the total variance of marginal treatment effects using the
delta method.

## Usage

``` r
marginal_variance(ald, ref_trt = NA, comp_trt = NA, scale, family)
```

## Arguments

- ald:

  Aggregate-level data

- ref_trt:

  Treatment labels reference (common; e.g. placebo)

- comp_trt:

  Treatment labels comparator

- scale:

  A scaling parameter for the calculation.

- family:

  A character string specifying the family distribution (e.g.,
  "binomial").

## Value

Numeric total variance of marginal treatment effects.

## Examples

``` r
ald <- data.frame(trt = c("B","C","B","C"),
                  variable = c(NA, NA, "y", "y"),
                  statistic = c("N", "N", "sum", "sum"),
                  value = c(100, 100, 50, 60))
                  
marginal_variance(ald, ref_trt = "C", comp_trt = "B",
                  scale = "log_odds", family = "binomial")
#> [1] 0.08166667
```

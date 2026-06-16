# Marginal treatment effect from reported event counts

Computes the relative treatment effect from aggregate-level data using
event counts.

## Usage

``` r
marginal_treatment_effect(ald, ref_trt = NA, comp_trt = NA, scale, family)
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

Numeric relative treatment effect.

## Examples

``` r
ald <- data.frame(trt = c("B","C","B","C"),
                  variable = c(NA, NA, "y", "y"),
                  statistic = c("N", "N", "sum", "sum"),
                  value = c(100, 100, 50, 60))
                  
marginal_treatment_effect(ald, ref_trt = "C", comp_trt = "B",
                          scale = "log_odds", family = "binomial")
#> [1] -0.4054651
```

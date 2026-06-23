<div id="main" class="col-md-9" role="main">

# Marginal treatment effect from reported event counts

<div class="ref-description section level2">

Computes the relative treatment effect from aggregate-level data using
event counts.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
marginal_treatment_effect(ald, ref_trt = NA, comp_trt = NA, scale, family)
```

</div>

</div>

<div class="section level2">

## Arguments

-   ald:

    Aggregate-level data

-   ref_trt:

    Treatment labels reference (common; e.g. placebo)

-   comp_trt:

    Treatment labels comparator

-   scale:

    A scaling parameter for the calculation.

-   family:

    A character string specifying the family distribution (e.g.,
    "binomial").

</div>

<div class="section level2">

## Value

Numeric relative treatment effect.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
ald <- data.frame(trt = c("B","C","B","C"),
                  variable = c(NA, NA, "y", "y"),
                  statistic = c("N", "N", "sum", "sum"),
                  value = c(100, 100, 50, 60))
                  
marginal_treatment_effect(ald, ref_trt = "C", comp_trt = "B",
                          scale = "log_odds", family = "binomial")
#> [1] -0.4054651
```

</div>

</div>

</div>

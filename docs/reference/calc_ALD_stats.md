<div id="main" class="col-md-9" role="main">

# Aggregate-level data mean and variance statistics

<div class="ref-description section level2">

Computes the mean and variance of marginal treatment effects for
aggregate-level trial data.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
calc_ALD_stats(strategy, analysis_params)
```

</div>

</div>

<div class="section level2">

## Arguments

-   strategy:

    A list containing the strategy details, including the family
    distribution.

-   analysis_params:

    A list containing:

    -   `ald` Aggregate-level trial data

    -   `ref_trt` Treatment labels reference (common; e.g. placebo)

    -   `ald_comp` Treatment labels comparator

    -   `scale` A scaling parameter for the calculation. From
        "log_odds", "risk_difference", "log_relative_risk".

</div>

<div class="section level2">

## Value

A list containing:

-   `mean`:

    The marginal treatment effect mean.

-   `var`:

    The marginal treatment effect variance.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`marginal_treatment_effect()` `marginal_variance()`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

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

</div>

</div>

</div>

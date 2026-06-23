<div id="main" class="col-md-9" role="main">

# Calculate trial variance

<div class="ref-description section level2">

Computes the variance of treatment effects for a trial based on the
specified family distribution.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
calculate_trial_variance(ald, tid, effect, family)
```

</div>

</div>

<div class="section level2">

## Arguments

-   ald:

    Aggregate-level data. Data frame.

-   tid:

    Treatment identifier used to extract relevant columns from `ald`.

-   effect:

    A character string specifying the effect scale (e.g., "log_odds",
    "risk_difference").

-   family:

    A character string specifying the model family (e.g., "binomial",
    "gaussian").

</div>

<div class="section level2">

## Value

Numeric computed variance of treatment effects.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
ald <- data.frame(trt = c("B","C","B","C"),
                  variable = c(NA, NA, "y", "y"),
                  statistic = c("N", "N", "sum", "sum"),
                  value = c(100, 100, 50, 60))
                  
calculate_trial_variance(
  ald, tid = "B", effect = "log_odds", family = "binomial")
#> [1] 0.04
  
```

</div>

</div>

</div>

<div id="main" class="col-md-9" role="main">

# Calculate Average Treatment Effect

<div class="ref-description section level2">

Computes the average treatment effect (ATE) based on the specified
effect scale.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
calculate_ate(mean_comp, mean_ref, effect)
```

</div>

</div>

<div class="section level2">

## Arguments

-   mean_comp, mean_ref:

    Mean of the outcome for the comparator and reference / common

-   effect:

    A character string specifying the effect scale. Options are:

    "log_odds"

    :   Log-odds difference.

    "risk_difference"

    :   Risk difference.

    "delta_z"

    :   Probit scale difference (z-scores).

    "log_relative_risk_rare_events"

    :   Log relative risk for rare events.

    "log_relative_risk"

    :   Log relative risk.

</div>

<div class="section level2">

## Value

Numeric computed average treatment effect on the specified scale.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
calculate_ate(mean_comp = 0.7, mean_ref = 0.5, effect = "log_odds")
#> [1] 0.8472979
calculate_ate(mean_comp = 0.7, mean_ref = 0.5, effect = "risk_difference")
#> [1] 0.2
```

</div>

</div>

</div>

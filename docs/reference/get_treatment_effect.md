<div id="main" class="col-md-9" role="main">

# Get treatment effect scale corresponding to a link function

<div class="ref-description section level2">

Maps a given link function to its corresponding treatment effect scale.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
get_treatment_effect(link)
```

</div>

</div>

<div class="section level2">

## Arguments

-   link:

    A character string specifying the link function. Options are:

    "logit"

    :   Log-odds scale.

    "identity"

    :   Risk difference.

    "probit"

    :   Probit scale.

    "cloglog"

    :   Log relative risk for rare events.

    "log"

    :   Log relative risk.

</div>

<div class="section level2">

## Value

A character string representing the treatment effect scale.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
 get_treatment_effect(link = "logit")
#> [1] "log_odds"
 get_treatment_effect(link = "identity")
#> [1] "mean_difference"
 
```

</div>

</div>

</div>

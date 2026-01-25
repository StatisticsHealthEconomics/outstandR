# Get treatment effect scale corresponding to a link function

Maps a given link function to its corresponding treatment effect scale.

## Usage

``` r
get_treatment_effect(link)
```

## Arguments

- link:

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

## Value

A character string representing the treatment effect scale.

## Examples

``` r
 get_treatment_effect(link = "logit")
#> [1] "log_odds"
 get_treatment_effect(link = "identity")
#> [1] "mean_difference"
 
```

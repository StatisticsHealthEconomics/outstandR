<div id="main" class="col-md-9" role="main">

# Variance estimate by pooling

<div class="ref-description section level2">

Use combining rules to estimate.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
var_by_pooling(n_imp, bar.v, b)
```

</div>

</div>

<div class="section level2">

## Arguments

-   n_imp:

    Number of syntheses used in analysis stage (high for low Monte Carlo
    error)

-   bar.v:

    "within" variance (average of variance point estimates)

-   b:

    "between" variance (sample variance of point estimates)

</div>

<div class="section level2">

## Value

Numeric value of variance estimate using pooling.

</div>

</div>

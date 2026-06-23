<div id="main" class="col-md-9" role="main">

# Wald-type interval estimates

<div class="ref-description section level2">

Constructed using t-distribution with nu degrees of freedom.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
wald_type_interval(n_imp, bar.v, b)
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

Numeric value of Wald-type interval estimates.

</div>

</div>

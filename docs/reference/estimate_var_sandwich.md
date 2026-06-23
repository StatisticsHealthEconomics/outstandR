<div id="main" class="col-md-9" role="main">

# Estimate Variance Sandwich Estimator

<div class="ref-description section level2">

Computes the robust (sandwich) variance estimator for the treatment
effect.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
estimate_var_sandwich(strategy, analysis_params, ...)
```

</div>

</div>

<div class="section level2">

## Arguments

-   strategy:

    An object of class `strategy` created by functions such as
    `strategy_maic()`, `strategy_stc()`, or `strategy_mim()`. Contains
    modelling details like the formula and family.

-   analysis_params:

    List of analysis parameters (ipd, ald, etc.)

-   ...:

    Additional arguments

</div>

<div class="section level2">

## Value

Numeric variance estimate for the treatment contrast

</div>

</div>

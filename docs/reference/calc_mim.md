<div id="main" class="col-md-9" role="main">

# Multiple imputation marginalization (MIM)

<div class="ref-description section level2">

Multiple imputation marginalization (MIM)

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
calc_mim(strategy, analysis_params, ...)
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

    List of analysis parameters. Must contain `ipd` and `ald`.

-   ...:

    Additional argument to pass to Stan model

</div>

<div class="section level2">

## Value

A list containing:

-   `means`: A list containing named vectors of posterior means (one per
    synthesis `n_imp`):

    -   Comparator means.

    -   Reference means.

-   `model`: A list containing:

    -   `fit`: The first-stage `rstanarm::stan_glm()` object.

    -   `hats.v`: Vector of variance point estimates for each synthesis.

    -   `n_imp`: Number of posterior prediction draws (syntheses).

    -   `rho`, `N`, `stan_args`: Strategy and model parameters.

</div>

</div>

<div id="main" class="col-md-9" role="main">

# Calculate MAIC

<div class="ref-description section level2">

Calculate MAIC

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
calc_maic(strategy, analysis_params)
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

    List of analysis parameters. Must contain `ipd` (individual patient
    data) and `ald` (aggregated lead data).

</div>

<div class="section level2">

## Value

A list containing:

-   `means`: A list containing:

    -   `A`: Bootstrap estimates for comparator treatment group "A".

    -   `C`: Bootstrap estimates for reference treatment group "C".

-   `model`: A list containing model diagnostics derived from the
    original data:

    -   `weights`: Vector of calculated weights for the patients in
        `ipd`.

    -   `ESS`: The Effective Sample Size.

</div>

</div>

<div id="main" class="col-md-9" role="main">

# Calculate simulated treatment comparison statistics

<div class="ref-description section level2">

Calculate simulated treatment comparison statistics

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
calc_stc(strategy, analysis_params, ...)
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
    data).

-   ...:

    Additional arguments.

</div>

<div class="section level2">

## Value

A list containing:

-   `means`: A list containing:

    -   `A`: Mean for comparator treatment group "A".

    -   `C`: Mean for reference treatment group "C".

-   `model`: The fitted `stats::glm()` object.

</div>

</div>

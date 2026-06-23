<div id="main" class="col-md-9" role="main">

# Get covariate names

<div class="ref-description section level2">

Get covariate names

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
get_covariate_names(formula)
```

</div>

</div>

<div class="section level2">

## Arguments

-   formula:

    Linear regression `formula` object. Prognostic factors (PF) are main
    effects and effect modifiers (EM) are interactions with the
    treatment variable, e.g., y \~ X1 + trt + trt:X2. For covariates as
    both PF and EM use `*` syntax.

</div>

<div class="section level2">

## Value

Covariate names character vector

</div>

</div>

<div id="main" class="col-md-9" role="main">

# Get effect modifiers

<div class="ref-description section level2">

Get effect modifiers

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
get_eff_mod_names(formula, trt_var = "trt")
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

-   trt_var:

    Treatment variable name; Default 'trt'.

</div>

<div class="section level2">

## Value

Effect modifiers strings names

</div>

</div>

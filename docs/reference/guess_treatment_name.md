<div id="main" class="col-md-9" role="main">

# Guess treatment name

<div class="ref-description section level2">

Does a variable appear more than once in interactions? If not then pick
first LHS interaction term. Finally, if there are no interactions then
pick last main effect term.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
guess_treatment_name(formula)
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

Treatment name string.

</div>

</div>

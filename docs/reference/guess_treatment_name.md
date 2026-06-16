# Guess treatment name

Does a variable appear more than once in interactions? If not then pick
first LHS interaction term. Finally, if there are no interactions then
pick last main effect term.

## Usage

``` r
guess_treatment_name(formula)
```

## Arguments

- formula:

  Linear regression `formula` object. Prognostic factors (PF) are main
  effects and effect modifiers (EM) are interactions with the treatment
  variable, e.g., y ~ X1 + trt + trt:X2. For covariates as both PF and
  EM use `*` syntax.

## Value

Treatment name string.

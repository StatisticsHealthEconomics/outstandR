# Get covariate names

Get covariate names

## Usage

``` r
get_covariate_names(formula)
```

## Arguments

- formula:

  Linear regression `formula` object. Prognostic factors (PF) are main
  effects and effect modifiers (EM) are interactions with the treatment
  variable, e.g., y ~ X1 + trt + trt:X2. For covariates as both PF and
  EM use `*` syntax.

## Value

Covariate names character vector

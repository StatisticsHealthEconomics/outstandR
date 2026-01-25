# Get effect modifiers

Get effect modifiers

## Usage

``` r
get_eff_mod_names(formula, trt_var = "trt")
```

## Arguments

- formula:

  Linear regression `formula` object. Prognostic factors (PF) are main
  effects and effect modifiers (EM) are interactions with the treatment
  variable, e.g., y ~ X1 + trt + trt:X2. For covariates as both PF and
  EM use `*` syntax.

- trt_var:

  Treatment variable name; Default 'trt'.

## Value

Effect modifiers strings names

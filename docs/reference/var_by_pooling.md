# Variance estimate by pooling

Use combining rules to estimate.

## Usage

``` r
var_by_pooling(n_imp, bar.v, b)
```

## Arguments

- n_imp:

  Number of syntheses used in analysis stage (high for low Monte Carlo
  error)

- bar.v:

  "within" variance (average of variance point estimates)

- b:

  "between" variance (sample variance of point estimates)

## Value

Numeric value of variance estimate using pooling.

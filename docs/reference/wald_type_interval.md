# Wald-type interval estimates

Constructed using t-distribution with nu degrees of freedom.

## Usage

``` r
wald_type_interval(n_imp, bar.v, b)
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

Numeric value of Wald-type interval estimates.

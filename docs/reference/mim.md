# Multiple imputation marginalization (MIM)

Multiple imputation marginalization (MIM)

## Usage

``` r
mim(formula, ipd, ald, M = 1000, n.chains = 2, warmup = 1000, iters = 4000)
```

## Arguments

- ipd:

  Individual-level data

- ald:

  Aggregate-level data

- M:

  Number of syntheses used in analysis stage (high for low Monte Carlo
  error)

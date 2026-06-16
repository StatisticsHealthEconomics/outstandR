# G-computation using Stan

Calculate draws of binary responses from posterior predictive
distribution from the Bayesian G-computation method using Hamiltonian
Monte Carlo.

## Usage

``` r
gcomp_bayes(formula = NULL, ipd, ald)
```

## Arguments

- formula:

  Linear regression `formula` object

- ipd:

  Individual-level data

- ald:

  Aggregate-level data

## Value

A list of \\y^\*\_A\\ and \\y^\*\_C\\ posterior predictions

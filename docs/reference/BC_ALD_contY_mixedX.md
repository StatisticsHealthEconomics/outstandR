# Aggregate level patient data for continuous outcome, mixed covariates

This data set contains summaries of simulated patient covariate and
outcome values. Corresponds to IPD data set.

## Usage

``` r
data(BC_ALD_contY_mixedX)
```

## Format

### `y ~ X1 + X3 + X4 + trt + trt:(X2 + X3 + X4)`

- variable:

  String covariate or outcome name. From X1, X2, X3, X4, y.

- statistic:

  String summary statistic name. From mean, sd, prob, sum, N

- value:

  Numeric value

- trt:

  Treatment (arm) name. From B, C

## Source

Simulated data

## References

Remiro‐Azocar A, Heath A, Baio G (2022)

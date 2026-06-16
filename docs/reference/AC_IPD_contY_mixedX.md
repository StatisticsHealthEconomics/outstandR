# Individual-level patient data for continuous outcome, mixed covariates

This data set contains simulated patient covariate and outcome values.
Corresponds to ALD data set.

## Usage

``` r
data(AC_IPD_contY_mixedX)
```

## Format

### `y ~ X1 + X3 + X4 + trt + trt:(X2 + X3 + X4)`

- id:

  Numeric unique identifier

- X1:

  Numeric prognostic factor continuous covariate

- X2:

  Numeric prognostic factor and effect modifier binary covariate

- X3:

  Numeric prognostic factor and effect modifier continuous covariate

- X4:

  Numeric effect modifier binary covariate

- trt:

  Factor treatment identifier. Levels A, C

- y:

  Integer binary outcome

- true_eta:

  Numeric linear predictor

## Source

Simulated data

## References

Remiro‐Azocar A, Heath A, Baio G (2022)

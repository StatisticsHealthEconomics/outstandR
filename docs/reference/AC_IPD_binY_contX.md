# Individual-level patient data for binary outcome, continuous covariates

This data set contains simulated patient covariate and outcome values.

## Usage

``` r
data(AC_IPD_binY_contX)
```

## Format

### `y ~ PF_cont_1 + PF_cont_2 + trt + trt:(EM_cont_1 + EM_cont_2)`

- id:

  Numeric unique identifier

- PF_cont_1:

  Numeric prognostic factor continuous covariate

- PF_cont_2:

  Numeric prognostic factor continuous covariate

- EM_cont_1:

  Numeric effect modifier continuous covariate

- EM_cont_2:

  Numeric effect modifier continuous covariate

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

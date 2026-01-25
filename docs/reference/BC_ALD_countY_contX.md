# Aggregate level patient data for count outcome, continuous covariates

This data set contains summaries of simulated patient covariate and
outcome values. Corresponds to IPD data set.

## Usage

``` r
data(BC_ALD_countY_contX)
```

## Format

### `y ~ PF_cont_1 + PF_cont_2 + trt + trt:(EM_cont_1 + EM_cont_2)`

- variable:

  String covariate or outcome name. From EM_cont_1, EM_cont_2,
  PF_cont_1, PF_cont_2, y.

- statistic:

  String summary statistic name. From mean, sd, sum, N

- value:

  Numeric value

- trt:

  Treatment (arm) name. From B, C

## Source

Simulated data

## References

Remiro‐Azocar A, Heath A, Baio G (2022)

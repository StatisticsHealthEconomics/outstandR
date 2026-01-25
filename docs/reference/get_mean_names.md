# Get mean names

Get mean names

## Usage

``` r
get_mean_names(ald, keep_nms)
```

## Arguments

- ald:

  Aggregate-level data. Single row matrix with summary statistics for
  each covariate and treatment outcomes. The format is 'mean.*' and
  'sd.*' for covariates and 'y.*.sum', 'y.*.bar', 'y.\*.sd' for
  treatments B and C. We assume a common distribution for each treatment
  arm.

- keep_nms:

  Variable names character vector

## Value

Mean names vector

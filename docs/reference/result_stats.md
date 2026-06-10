# Calculate and arrange result statistics

Combining output from aggregate level data studies BC and adjusted
individual level data studies AC into a single object.

## Usage

``` r
result_stats(ipd_stats, ald_stats, CI = 0.95)
```

## Arguments

- ipd_stats, ald_stats:

  Trial data for IPD and ALD

- CI:

  Confidence interval level, i.e. 1-alpha; default 0.95

## Value

List of ITC output:

- `contrasts`: A list for relative effects containing:

  - `means`

  - `variances`

  - `CI`

- `absolute`: A list for absolute effects containing:

  - `means`

  - `variances`

  - `CI`

# Convert aggregate data from long to wide format

Convert aggregate data from long to wide format

## Usage

``` r
reshape_ald_to_wide(df)
```

## Arguments

- df:

  A Dataframe of ALD

## Value

Data frame in wide format

## Examples

``` r
df <- 
  data.frame(
    variable = c("age", "age", "y", "y", "y", "y", "y", "y", "y", "y"),
    statistic = c("mean", "sd", "sum", "bar", "sd", "N", "sum", "bar", "sd", "N"),
    trt = c(NA, NA, "B", "B", "B", "B", "C", "C", "C", "C"),
    value = c(1,1,1,1,1,1,1,1,1,1))
```

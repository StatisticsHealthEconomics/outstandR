<div id="main" class="col-md-9" role="main">

# Convert aggregate data from long to wide format

<div class="ref-description section level2">

Convert aggregate data from long to wide format

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
reshape_ald_to_wide(df)
```

</div>

</div>

<div class="section level2">

## Arguments

-   df:

    A Dataframe of ALD

</div>

<div class="section level2">

## Value

Data frame in wide format

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
df <- 
  data.frame(
    variable = c("age", "age", "y", "y", "y", "y", "y", "y", "y", "y"),
    statistic = c("mean", "sd", "sum", "bar", "sd", "N", "sum", "bar", "sd", "N"),
    trt = c(NA, NA, "B", "B", "B", "B", "C", "C", "C", "C"),
    value = c(1,1,1,1,1,1,1,1,1,1))
```

</div>

</div>

</div>

<div id="main" class="col-md-9" role="main">

# Aggregate level patient data for count outcome, continuous covariates

<div class="ref-description section level2">

This data set contains summaries of simulated patient covariate and
outcome values. Corresponds to IPD data set.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
data(BC_ALD_countY_contX)
```

</div>

</div>

<div class="section level2">

## Format

<div class="section">

### `y ~ PF_cont_1 + PF_cont_2 + trt + trt:(EM_cont_1 + EM_cont_2)`

-   variable:

    String covariate or outcome name. From EM_cont_1, EM_cont_2,
    PF_cont_1, PF_cont_2, y.

-   statistic:

    String summary statistic name. From mean, sd, sum, N

-   value:

    Numeric value

-   trt:

    Treatment (arm) name. From B, C

</div>

</div>

<div class="section level2">

## Source

Simulated data

</div>

<div class="section level2">

## References

Remiro‐Azocar A, Heath A, Baio G (2022)

</div>

</div>

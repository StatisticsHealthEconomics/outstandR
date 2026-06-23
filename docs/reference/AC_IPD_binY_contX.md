<div id="main" class="col-md-9" role="main">

# Individual-level patient data for binary outcome, continuous covariates

<div class="ref-description section level2">

This data set contains simulated patient covariate and outcome values.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
data(AC_IPD_binY_contX)
```

</div>

</div>

<div class="section level2">

## Format

<div class="section">

### `y ~ PF_cont_1 + PF_cont_2 + trt + trt:(EM_cont_1 + EM_cont_2)`

-   id:

    Numeric unique identifier

-   PF_cont_1:

    Numeric prognostic factor continuous covariate

-   PF_cont_2:

    Numeric prognostic factor continuous covariate

-   EM_cont_1:

    Numeric effect modifier continuous covariate

-   EM_cont_2:

    Numeric effect modifier continuous covariate

-   trt:

    Factor treatment identifier. Levels A, C

-   y:

    Integer binary outcome

-   true_eta:

    Numeric linear predictor

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

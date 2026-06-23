<div id="main" class="col-md-9" role="main">

# Individual-level patient data for continuous outcome, mixed covariates

<div class="ref-description section level2">

This data set contains simulated patient covariate and outcome values.
Corresponds to ALD data set.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
data(AC_IPD_contY_mixedX)
```

</div>

</div>

<div class="section level2">

## Format

<div class="section">

### `y ~ X1 + X3 + X4 + trt + trt:(X2 + X3 + X4)`

-   id:

    Numeric unique identifier

-   X1:

    Numeric prognostic factor continuous covariate

-   X2:

    Numeric prognostic factor and effect modifier binary covariate

-   X3:

    Numeric prognostic factor and effect modifier continuous covariate

-   X4:

    Numeric effect modifier binary covariate

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

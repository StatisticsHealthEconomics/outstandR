<div id="main" class="col-md-9" role="main">

# MAIC bootstrap sample

<div class="ref-description section level2">

Matching-adjusted indirect comparison bootstrap sampling.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
maic.boot(
  data,
  indices,
  balance_matrix,
  outcome_x_matrix,
  outcome_y,
  ald_targets,
  scaling_factors,
  trt_var,
  family,
  hat_w = NULL,
  ipd = NULL,
  outcome_model = NULL,
  balance_model = NULL,
  ald = NULL,
  moments = 1,
  int = FALSE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   data:

    Individual-level patient data (data frame).

-   indices:

    Vector of indices, same length as original, which define the
    bootstrap sample.

-   balance_matrix:

    Pre-computed balance matrix.

-   outcome_x\_matrix:

    Pre-computed outcome design matrix.

-   outcome_y:

    Pre-computed outcome vector.

-   ald_targets:

    Vector of ALD targets.

-   scaling_factors:

    Vector of scaling factors.

-   trt_var:

    Treatment variable name.

-   family:

    A 'family' object specifying the distribution and link function.

-   hat_w:

    MAIC weights; default `NULL` which calls `maic_weights()`.

-   ipd:

    Backwards compatibility IPD data (optional).

-   outcome_model:

    Backwards compatibility outcome model formula (optional).

-   balance_model:

    Backwards compatibility balance model formula (optional).

-   ald:

    Backwards compatibility ALD data (optional).

-   moments:

    Backwards compatibility moments (default 1).

-   int:

    Backwards compatibility interactions flag (default FALSE).

</div>

<div class="section level2">

## Value

Vector of fitted probabilities for treatments *A* and *C*

</div>

<div class="section level2">

## See also

<div class="dont-index">

`calc_IPD_stats.maic()`

</div>

</div>

</div>

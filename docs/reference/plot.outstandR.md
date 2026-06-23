<div id="main" class="col-md-9" role="main">

# Default Plot Method for outstandR Objects

<div class="ref-description section level2">

Default Plot Method for outstandR Objects

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
# S3 method for class 'outstandR'
plot(
  x,
  ...,
  type = c("both", "contrasts", "absolute"),
  labels = NULL,
  include_naive = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    An object of class 'outstandR' or a list of 'outstandR' objects.

-   ...:

    Additional 'outstandR' objects for comparison.

-   type:

    Character, one of "both" (default), "contrasts", or "absolute".

-   labels:

    Optional character vector of names for the models.

-   include_naive:

    Logical. Should naive (unadjusted) estimates be included in the
    plot? Default is TRUE.

</div>

<div class="section level2">

## Value

A `ggplot2::ggplot()` object representing the forest plot of the
results.

</div>

</div>

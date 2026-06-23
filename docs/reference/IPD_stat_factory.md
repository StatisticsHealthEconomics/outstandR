<div id="main" class="col-md-9" role="main">

# Factory function for creating calc_IPD_stats methods

<div class="ref-description section level2">

Creates a method for computing IPD mean and variance statistics based on
the supplied function.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
IPD_stat_factory(ipd_fun)
```

</div>

</div>

<div class="section level2">

## Arguments

-   ipd_fun:

    A function that computes mean and variance statistics for
    individual-level patient data.

</div>

<div class="section level2">

## Value

A function that computes mean and variance statistics for a given
strategy.

</div>

</div>

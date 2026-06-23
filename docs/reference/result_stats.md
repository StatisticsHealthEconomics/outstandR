<div id="main" class="col-md-9" role="main">

# Calculate and arrange result statistics

<div class="ref-description section level2">

Combining output from aggregate level data studies BC and adjusted
individual level data studies AC into a single object.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
result_stats(ipd_stats, ald_stats, CI = 0.95)
```

</div>

</div>

<div class="section level2">

## Arguments

-   ipd_stats, ald_stats:

    Trial data for IPD and ALD

-   CI:

    Confidence interval level, i.e. 1-alpha; default 0.95

</div>

<div class="section level2">

## Value

List of ITC output:

-   `contrasts`: A list for relative effects containing:

    -   `means`

    -   `variances`

    -   `CI`

-   `absolute`: A list for absolute effects containing:

    -   `means`

    -   `variances`

    -   `CI`

</div>

</div>

<div id="main" class="col-md-9" role="main">

# Calculate Naive (Unadjusted) Bucher Indirect Comparison

<div class="ref-description section level2">

Calculates the unadjusted indirect comparison between two active
treatments via a common comparator using the standard Bucher method.
Estimates are calculated on the link scale (e.g., log-odds) and
back-transformed to the requested scale.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
calc_bucher_naive(
  ipd_trial,
  ald_trial,
  outcome_model,
  family,
  ref_trt = "C",
  scale = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   ipd_trial:

    Data frame containing the Individual Patient Data.

-   ald_trial:

    Data frame containing the Aggregate Level Data.

-   outcome_model:

    Formula specifying the outcome model (used to identify variables).

-   family:

    A family object (e.g., `binomial(link="logit")`).

-   ref_trt:

    Character string identifying the common comparator treatment.

-   scale:

    Character string specifying the desired output scale (e.g.,
    "risk_difference").

</div>

<div class="section level2">

## Value

A list containing the naive point estimate, variance, and standard
error.

</div>

</div>

<div id="main" class="col-md-9" role="main">

# Estimate MAIC weights

<div class="ref-description section level2">

Matching-adjusted indirect comparison weights. Method is taken from
(Signorovitch et al. 2010) .

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
maic_weights(X_EM)
```

</div>

</div>

<div class="section level2">

## Arguments

-   X_EM:

    Centred \\(S=1\\) effect modifiers IPD covariates; matrix or data
    frame

</div>

<div class="section level2">

## Value

Estimated weights for each individual; vector

</div>

<div class="section level2">

## References

Signorovitch JE, Wu EQ, Yu AP, Gerrits CM, Kantor E, Bao Y, Gupta SR,
Mulani PM (2010). “Comparative effectiveness without head-to-head
trials: A method for matching-adjusted indirect comparisons applied to
psoriasis treatment with adalimumab or etanercept.” *Pharmacoeconomics*,
**28**(10), 935–945. ISSN 11707690,
[doi:10.2165/11538370-000000000-00000](https://doi.org/10.2165/11538370-000000000-00000)
, <https://pubmed.ncbi.nlm.nih.gov/20831302/>.

</div>

</div>

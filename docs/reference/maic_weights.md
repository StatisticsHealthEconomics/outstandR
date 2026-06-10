# Estimate MAIC weights

Matching-adjusted indirect comparison weights. Method is taken from
(Signorovitch et al. 2010) .

## Usage

``` r
maic_weights(X_EM)
```

## Arguments

- X_EM:

  Centred \\S=1\\ effect modifiers IPD covariates; matrix or data frame

## Value

Estimated weights for each individual; vector

## References

Signorovitch JE, Wu EQ, Yu AP, Gerrits CM, Kantor E, Bao Y, Gupta SR,
Mulani PM (2010). “Comparative effectiveness without head-to-head
trials: A method for matching-adjusted indirect comparisons applied to
psoriasis treatment with adalimumab or etanercept.” *Pharmacoeconomics*,
**28**(10), 935–945. ISSN 11707690.
[doi:10.2165/11538370-000000000-00000](https://doi.org/10.2165/11538370-000000000-00000)
. <https://pubmed.ncbi.nlm.nih.gov/20831302/>.

# Compute Correlations for Multivariate Mom Model

Compute Correlations for Multivariate Mom Model

## Usage

``` r
compute_mom_multivariate_correlations(
  data,
  responses,
  random_facets,
  aov_models,
  aov_formulas,
  random_effect_cors = list(),
  compute_rescor = FALSE,
  compute_re_cors = FALSE
)
```

## Arguments

- data:

  The data frame

- responses:

  Character vector of response variable names

- random_facets:

  Character vector of random effect facet names

- aov_models:

  List of fitted aov models (one per response variable)

- random_effect_cors:

  List of random effect correlation specifications

- compute_rescor:

  Whether to compute residual correlations

- compute_re_cors:

  Whether to compute random effect correlations

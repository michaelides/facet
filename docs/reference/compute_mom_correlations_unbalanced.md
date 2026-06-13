# Compute Correlations for Unbalanced Multivariate Mom Model

Computes residual and random effect correlations using pairwise complete
cases for unbalanced multivariate designs.

## Usage

``` r
compute_mom_correlations_unbalanced(
  data,
  responses,
  random_facets,
  aov_models,
  data_per_dim,
  aov_formulas = list(),
  random_effect_cors = list(),
  compute_rescor = FALSE,
  compute_re_cors = FALSE
)
```

## Arguments

- data:

  Original data frame

- responses:

  Character vector of response variable names

- random_facets:

  Character vector of random effect facet names

- aov_models:

  List of fitted aov models (one per response variable)

- data_per_dim:

  List of data frames used for each dimension

- aov_formulas:

  List of aov formula strings per response (with Error terms)

- random_effect_cors:

  List of random effect correlation specifications

- compute_rescor:

  Whether to compute residual correlations

- compute_re_cors:

  Whether to compute random effect correlations

## Value

List with residual and random effect correlations

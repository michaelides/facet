# Extract Covariances from brms Multivariate Model

Extracts residual covariances and random effect covariances from a
multivariate brms model by computing them from the full posterior.

## Usage

``` r
extract_covariances_brms(model)
```

## Arguments

- model:

  A brmsfit object from a multivariate model.

## Value

A list with:

- residual_cov:

  Tibble with residual covariances (dim1, dim2, estimate, se, lower,
  upper, Rhat, Bulk_ESS, Tail_ESS)

- random_effect_cov:

  Named list of tibbles for each facet with correlated random effects

- residual_cov_matrix:

  Matrix of residual covariance point estimates

- random_effect_cov_matrix:

  Named list of matrices for random effect covariances

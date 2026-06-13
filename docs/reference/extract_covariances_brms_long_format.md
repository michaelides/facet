# Extract Covariances from brms Long-Format Multivariate Model

Extracts residual covariances and random effect covariances from a
long-format multivariate brms model. Converts correlations to
covariances using Cov = Cor \* SD_i \* SD_j.

## Usage

``` r
extract_covariances_brms_long_format(model, dimension_var, data, vc = NULL)
```

## Arguments

- model:

  A fitted brms model

- dimension_var:

  Name of the dimension variable in the data

- data:

  The original data used to fit the model

- vc:

  Optional variance components tibble (for extracting SDs)

## Value

A list with:

- residual_cov: Tibble with residual covariances (dim1, dim2, estimate,
  se, lower, upper)

- random_effect_cov: Named list of tibbles for each facet with
  covariances

- residual_cov_matrix: Matrix of residual covariance point estimates

- random_effect_cov_matrix: Named list of matrices for random effect
  covariances

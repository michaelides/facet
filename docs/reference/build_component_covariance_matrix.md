# Build Variance-Covariance Matrix for a Component

Constructs a variance-covariance matrix for a specific variance
component in a multivariate design. The diagonal elements are the
variances for each dimension, and off-diagonal elements are covariances
(when available).

## Usage

``` r
build_component_covariance_matrix(
  vc,
  component,
  dimensions,
  correlations = NULL,
  scale_factor = 1
)
```

## Arguments

- vc:

  A variance components tibble with columns: component, dim, var

- component:

  Character string naming the component (e.g., "Person", "Residual")

- dimensions:

  Character vector of dimension names

- correlations:

  List containing covariance information from the G-study:

  - `residual_cov`: Tibble of residual covariances

  - `random_effect_cov`: Named list of random effect covariances per
    facet

## Value

A symmetric variance-covariance matrix with dimensions as row/column
names

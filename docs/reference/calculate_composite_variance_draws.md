# Calculate Composite Variance Draws

Computes posterior draws of composite variance for each component type
using weighted variance-covariance matrices.

## Usage

``` r
calculate_composite_variance_draws(
  vc_draws,
  cov_draws,
  weights,
  scale_factor = 1
)
```

## Arguments

- vc_draws:

  Named list (by dimension) of named lists (by component) of variance
  draws

- cov_draws:

  Named list from extract_covariance_draws()

- weights:

  Named numeric vector of weights (names = dimensions)

- scale_factor:

  Factor to scale variances and covariances (for D-study scaling)

## Value

Named list (by component) of posterior draws for composite variance

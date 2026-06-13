# Compute Composite Error Variance

Calculates the composite relative or absolute error variance by summing
weighted variance-covariance matrices for all error components.

## Usage

``` r
compute_composite_error_variance(
  vc,
  dimensions,
  weights,
  correlations,
  universe_spec,
  error_spec,
  aggregation,
  residual_is,
  error_type,
  n = NULL,
  object_spec = NULL
)
```

## Arguments

- vc:

  A variance components tibble

- dimensions:

  Character vector of dimension names

- weights:

  Named numeric vector of weights

- correlations:

  List of covariance information from the G-study

- universe_spec:

  Character vector of universe components

- error_spec:

  Optional character vector of error components

- aggregation:

  Optional aggregation specification

- residual_is:

  Optional residual specification

- error_type:

  Either "relative" or "absolute"

- scale_factor:

  Factor by which to scale covariances

## Value

The composite error variance as a numeric scalar

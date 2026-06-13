# Compute Absolute Error Variance from Draws

Compute Absolute Error Variance from Draws

## Usage

``` r
compute_absolute_error_draws(
  scaled_draws,
  universe_spec,
  error_spec,
  agg_facets,
  residual_is = NULL
)
```

## Arguments

- scaled_draws:

  Named list of scaled variance draws.

- universe_spec:

  Character vector of universe components.

- error_spec:

  Character vector of error components (or NULL).

- agg_facets:

  Character vector of aggregation facets.

## Value

Numeric vector of absolute error variances.

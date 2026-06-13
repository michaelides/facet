# Compute Relative Error Variance from Draws

Compute Relative Error Variance from Draws

## Usage

``` r
compute_relative_error_draws(
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

Numeric vector of relative error variances.

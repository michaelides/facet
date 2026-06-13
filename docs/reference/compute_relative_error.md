# Compute Relative Error Variance

Calculates the relative error variance for a D-study design. Relative
error includes variance components that affect relative rankings
(interactions with the universe components).

## Usage

``` r
compute_relative_error(
  vc,
  universe_spec,
  error_spec = NULL,
  agg_facets = NULL,
  residual_is = NULL
)
```

## Arguments

- vc:

  Variance components tibble.

- universe_spec:

  Character vector of universe component names.

- error_spec:

  Character vector of error component names, or NULL for default.

## Value

The relative error variance.

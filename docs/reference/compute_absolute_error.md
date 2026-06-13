# Compute Absolute Error Variance

Calculates the absolute error variance for a D-study design. Absolute
error includes all variance components that affect absolute scores (all
components except the universe components).

## Usage

``` r
compute_absolute_error(
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

The absolute error variance.

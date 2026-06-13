# Identify Error Components for D-Study (Internal)

Determines which variance components belong to relative and absolute
error. This consolidates logic used by both point-estimate and
draw-based computations.

## Usage

``` r
identify_error_components_for_draws(
  components,
  universe_spec,
  error_spec = NULL,
  agg_facets = NULL
)
```

## Arguments

- components:

  Character vector of component names

- universe_spec:

  Character vector of universe components

- error_spec:

  Character vector of explicit error components (or NULL)

- agg_facets:

  Character vector of aggregation facets (or NULL)

## Value

List with `is_relative_error` and `is_absolute_error` logical vectors

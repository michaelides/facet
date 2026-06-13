# Classify Error Components for D-Study

Unified function that determines which variance components belong to
relative and absolute error, consolidating the logic previously
duplicated across
[`compute_relative_error()`](https://github.com/yourorg/facet/reference/compute_relative_error.md),
[`compute_absolute_error()`](https://github.com/yourorg/facet/reference/compute_absolute_error.md),
and
[`identify_error_components_for_draws()`](https://github.com/yourorg/facet/reference/identify_error_components_for_draws.md).

## Usage

``` r
classify_error_components(
  components,
  universe_spec,
  error_spec = NULL,
  agg_facets = NULL
)
```

## Arguments

- components:

  Character vector of component names.

- universe_spec:

  Character vector of universe components.

- error_spec:

  Character vector of explicit error components (or NULL for default).

- agg_facets:

  Character vector of aggregation facets (or NULL).

## Value

A list with:

- is_relative_error:

  Logical vector: TRUE for components in relative error

- is_absolute_error:

  Logical vector: TRUE for components in absolute error

- error_spec_used:

  Logical: TRUE if an explicit error_spec was provided

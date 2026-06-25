# Compute Scale Factor for D-Study Variance Component

Calculates the scaling factor for a variance component based on the
D-study sample sizes. For main effects, divides by n_facet. For
interactions, divides by the product of n for each non-object facet.
Facets in `object_spec` are excluded from the divisor (the object of
measurement is not averaged over).

## Usage

``` r
compute_scale_factor_from_facets(facets, n, object_spec = NULL)
```

## Arguments

- facets:

  Character vector of facet names in the component.

- n:

  Named list of sample sizes for each facet.

- object_spec:

  Character vector of object component names. Facets in `object_spec`
  are excluded from the divisor. Default NULL (no exclusion).

## Value

The scale factor (numeric).

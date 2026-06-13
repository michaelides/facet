# Compute Scale Factor for D-Study Variance Component

Calculates the scaling factor for a variance component based on the
D-study sample sizes. For main effects, divides by n_facet. For
interactions, divides by the product of n for each facet.

## Usage

``` r
compute_scale_factor_from_facets(facets, n)
```

## Arguments

- facets:

  Character vector of facet names in the component.

- n:

  Named list of sample sizes for each facet.

## Value

The scale factor (numeric).

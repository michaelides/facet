# Reorder Variance Components to Match Formula Order

Reorders the variance components tibble to follow the user-specified
formula order, with Residual always last. Also renames interaction
components to match the user's specification (e.g., "Rater:Task" vs
"Task:Rater").

## Usage

``` r
reorder_variance_components(vc, facet_specs)
```

## Arguments

- vc:

  A tibble of variance components

- facet_specs:

  Character vector of facet specifications in user order

## Value

The vc tibble reordered to match facet_specs, with Residual last

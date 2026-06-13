# Parse Facet Names from Component String

Extracts facet names from a variance component string. Handles
interactions (e.g., "p:r:s"), main effects, and residual.

## Usage

``` r
parse_component_facets(component)
```

## Arguments

- component:

  Character string naming the variance component.

## Value

A character vector of facet names.

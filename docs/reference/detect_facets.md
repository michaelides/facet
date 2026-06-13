# Detect Facets from a Formula

Identifies which variables in the formula are facets (sources of
variance) in the generalizability theory sense.

## Usage

``` r
detect_facets(formula, data = NULL)
```

## Arguments

- formula:

  A formula object.

- data:

  A data frame (optional, used for validation).

## Value

A character vector of facet names.

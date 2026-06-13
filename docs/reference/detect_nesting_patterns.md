# Detect Nesting Patterns from Data

Analyzes the data structure to determine if facets are nested. A facet A
is considered nested in facet B if each level of A appears with only one
level of B.

## Usage

``` r
detect_nesting_patterns(data, facets)
```

## Arguments

- data:

  A data frame.

- facets:

  Character vector of facet names to analyze.

## Value

A list where each element describes a nesting relationship:

- nested_facet:

  Name of the nested facet

- nesting_facet:

  Name of the facet that the nested facet is within

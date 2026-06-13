# Parse Formula to Extract Facets

Extracts information about the response variable, fixed effects, and
random effects from a G-study formula.

## Usage

``` r
parse_g_formula(formula)
```

## Arguments

- formula:

  A formula object.

## Value

A list with components:

- response:

  Character string naming the response variable

- fixed:

  Character vector of fixed effect terms

- random:

  Character vector of random effect terms (as parsed)

- random_facets:

  Character vector of individual facet names from random effects

- random_facet_specs:

  Character vector of facet specifications as user specified (e.g.,
  "Rater:Task")

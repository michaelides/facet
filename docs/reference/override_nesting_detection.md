# Override Nesting Detection with User Specification

Merges auto-detected nesting with user-specified nesting relationships.

## Usage

``` r
override_nesting_detection(detected, user_nested, all_facets)
```

## Arguments

- detected:

  List of auto-detected nesting relationships.

- user_nested:

  Named list of user-specified nesting (nested = nesting).

- all_facets:

  Character vector of all facet names.

## Value

Updated list of nesting relationships.

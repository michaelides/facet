# Validate D-Study Inputs and Extract Basic Properties

Validates the gstudy object, determines if multivariate, and processes
weights.

## Usage

``` r
validate_dstudy_inputs(gstudy_obj, weights = NULL)
```

## Arguments

- gstudy_obj:

  A gstudy or mgstudy object

- weights:

  Optional numeric vector of weights

## Value

A list with is_multivariate, dimensions, weights

# Optimize Weights for Subscale VAR

Main function for "subscale" optimization method.

## Usage

``` r
optimize_subscale_weights_internal(
  dstudy_obj,
  subscale = NULL,
  optimize_target = "rel"
)
```

## Arguments

- dstudy_obj:

  A dstudy object

- subscale:

  Target subscale name (NULL for all)

- optimize_target:

  "rel" or "abs"

## Value

List with optimization results

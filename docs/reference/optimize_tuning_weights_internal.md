# Optimize Weights via Grid Search (Tuning)

Main function for "tuning" optimization method.

## Usage

``` r
optimize_tuning_weights_internal(
  dstudy_obj,
  grid_resolution = 0.1,
  optimize_target = "rel"
)
```

## Arguments

- dstudy_obj:

  A dstudy object

- grid_resolution:

  Step size for grid

- optimize_target:

  "rel" or "abs"

## Value

List with optimization results

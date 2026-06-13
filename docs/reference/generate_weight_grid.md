# Generate Weight Grid for Tuning Method

Creates a grid of weight combinations for exhaustive search.

## Usage

``` r
generate_weight_grid(k, resolution = 0.1, min_weight = 0.01)
```

## Arguments

- k:

  Number of dimensions

- resolution:

  Step size for grid (default 0.1)

- min_weight:

  Minimum weight allowed (default 0.01)

## Value

Tibble with all valid weight combinations

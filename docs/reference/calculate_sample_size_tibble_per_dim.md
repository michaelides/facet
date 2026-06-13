# Calculate Sample Size Tibble Per Dimension

Converts per-dimension sample size information into a tibble with a dim
column.

## Usage

``` r
calculate_sample_size_tibble_per_dim(sample_size_info_per_dim)
```

## Arguments

- sample_size_info_per_dim:

  List of sample_size_info objects per dimension.

## Value

A tibble with columns: dim, effect, type, n

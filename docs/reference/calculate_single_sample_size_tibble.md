# Convert Sample Size Info to Tibble Format

Converts the list-based sample_size_info into a compact tibble format
suitable for display in print and summary methods.

## Usage

``` r
calculate_single_sample_size_tibble(sample_size_info)
```

## Arguments

- sample_size_info:

  List from calculate_sample_size_info().

## Value

A tibble with columns:

- effect:

  Name of the effect (facet, interaction, or nested description)

- type:

  Type of effect: "main", "interaction", "residual", or "nested"

- n:

  Sample size for the effect

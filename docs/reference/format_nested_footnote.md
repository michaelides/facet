# Format Nested Details as Footnote

Creates formatted footnote text for nested effects, showing group
counts, means, and harmonic means when applicable.

## Usage

``` r
format_nested_footnote(nested_info)
```

## Arguments

- nested_info:

  List of nested effect information from
  calculate_nested_sample_sizes().

## Value

Character vector of formatted footnote lines, one per nested effect.
Returns character(0) if nested_info is NULL or empty.

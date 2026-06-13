# Summarize Correlation Tibble for Display

Formats a correlation tibble for printing, similar to summarize_vc.

## Usage

``` r
summarize_cor(cor_tibble, digits = 3)
```

## Arguments

- cor_tibble:

  A tibble with columns: dim1, dim2, estimate, se, lower, upper, Rhat,
  Bulk_ESS, Tail_ESS

- digits:

  Number of digits for rounding (default 3).

## Value

A tibble formatted for display.

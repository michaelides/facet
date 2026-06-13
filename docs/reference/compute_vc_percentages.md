# Compute Variance Component Percentages

Calculates the percentage of total variance accounted for by each
component.

## Usage

``` r
compute_vc_percentages(vc)
```

## Arguments

- vc:

  A tibble of variance components with a "var" column.

## Value

The input tibble with an added "pct" column.

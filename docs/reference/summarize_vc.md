# Summarize Variance Components

Creates a summary table of variance components suitable for reporting.

## Usage

``` r
summarize_vc(vc, digits = 4, scale = c("variance", "sd"))
```

## Arguments

- vc:

  A tibble of variance components.

- digits:

  Number of decimal places for rounding (default 4).

- scale:

  Scale for displaying results: "variance" (default) or "sd".

## Value

A formatted tibble.

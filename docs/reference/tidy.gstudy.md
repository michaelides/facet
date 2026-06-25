# Tidy Method for gstudy Objects

Returns a tidy tibble of variance components with numeric columns
rounded to `digits` decimal places (default 4).

Returns a tidy tibble of variance components for multivariate G-studies
with numeric columns rounded to `digits` decimal places (default 4).

## Usage

``` r
# S3 method for class 'gstudy'
tidy(x, digits = 4, ...)

# S3 method for class 'mgstudy'
tidy(x, digits = 4, ...)
```

## Arguments

- x:

  An mgstudy object.

- digits:

  Number of decimal places for rounding numeric columns (default 4).

- ...:

  Additional arguments (ignored).

## Value

A tibble.

A tibble.

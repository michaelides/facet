# Tidy Method for dstudy Objects

Returns a tidy tibble of D-study coefficients with numeric columns
rounded to `digits` decimal places (default 4).

## Usage

``` r
# S3 method for class 'dstudy'
tidy(x, digits = 4, ...)
```

## Arguments

- x:

  A dstudy object.

- digits:

  Number of decimal places for rounding numeric columns (default 4).

- ...:

  Additional arguments (ignored).

## Value

A tibble.

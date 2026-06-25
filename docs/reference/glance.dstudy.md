# Glance Method for dstudy Objects

Returns the unscaled variance component estimates from a D-study as a
tibble. For univariate D-studies the result is a one-row tibble with one
column per variance component, named `var_unscaled_<component>`. For
multivariate D-studies the result is a multi-row tibble with one row per
variance component and one column per dimension. Covariances are not
included; "Composite" rows from posterior estimation are dropped.
Numeric values are rounded to `digits` decimal places (default 4).

## Usage

``` r
# S3 method for class 'dstudy'
glance(x, digits = 4, ...)
```

## Arguments

- x:

  A dstudy object.

- digits:

  Number of decimal places for rounding numeric values (default 4).

- ...:

  Additional arguments (ignored).

## Value

A tibble of unscaled variance component estimates.

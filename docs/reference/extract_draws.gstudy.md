# Extract Draws from gstudy Objects

Extract posterior draws from a gstudy object that was fit with the brms
backend. This is a wrapper around brms::as_draws_matrix.

## Usage

``` r
extract_draws.gstudy(object, ...)
```

## Arguments

- object:

  A gstudy object.

- ...:

  Additional arguments passed to brms::as_draws_matrix.

## Value

A posterior draws_matrix object containing posterior draws for all model
parameters.

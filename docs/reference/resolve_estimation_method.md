# Resolve Estimation Method

Determines the appropriate estimation method based on the backend. Warns
and refits if posterior requested with non-brms backend. Warns and
overrides if simple requested with brms backend.

## Usage

``` r
resolve_estimation_method(gstudy_obj, estimation = NULL)
```

## Arguments

- gstudy_obj:

  A gstudy or mgstudy object

- estimation:

  Requested estimation method (NULL, "simple", or "posterior")

## Value

A list with gstudy_obj (possibly refit) and estimation (resolved)

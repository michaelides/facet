# Bind response variables for multivariate models

This function is a wrapper around
[`brms::mvbind()`](https://paulbuerkner.com/brms/reference/mvbind.html)
to allow specifying multivariate models via the brms estimator in
gstudy. See
[`brms::mvbind()`](https://paulbuerkner.com/brms/reference/mvbind.html)
for full documentation.

## Usage

``` r
mvbind(...)
```

## Arguments

- ...:

  Unquoted names of variables to bind.

## Value

A matrix-like object for brms formulas.

# Set residual correlation between response variables

This function is a wrapper around
[`brms::set_rescor()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)
to allow specifying residual correlations in multivariate models fit via
the brms backend. See
[`brms::set_rescor()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)
for full documentation.

## Usage

``` r
set_rescor(rescor)
```

## Arguments

- rescor:

  Logical; whether to estimate residual correlations.

## Value

A brmsformula helper object.

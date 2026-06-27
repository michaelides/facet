# Set up brms formulas

This function is a wrapper around
[`brms::bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
to allow specifying formulas for Bayesian models fit via the brms
estimator in gstudy. See
[`brms::bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
for full documentation.

## Usage

``` r
bf(formula, ...)
```

## Arguments

- formula:

  A formula object.

- ...:

  Additional arguments passed to
  [`brms::bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html).

## Value

A brmsformula object.

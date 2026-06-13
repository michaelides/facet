# Fit a Model Using brms

Fits a Bayesian mixed effects model using brms::brm(). Handles both
univariate and multivariate formulas.

## Usage

``` r
fit_brms(formula, data, prior = NULL, ...)
```

## Arguments

- formula:

  A formula for the model (can be brmsformula).

- data:

  A data frame.

- prior:

  A brmsprior object or list of priors as created by
  [`set_prior()`](https://github.com/yourorg/facet/reference/set_prior.md)
  or related functions. If NULL, default priors are used.

- ...:

  Additional arguments passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

## Value

A fitted brmsfit object.

# Extract Variance Component Draws from brms Model

Extracts posterior draws of variance components (SD^2) from a brms
model. For multivariate models, returns a nested list structure by
dimension.

## Usage

``` r
extract_variance_draws(gstudy_obj, draws)
```

## Arguments

- gstudy_obj:

  A gstudy object with brms estimator.

- draws:

  A posterior draws_matrix object (e.g., from brms::as_draws_matrix()).

## Value

A named list. For univariate: list of variance draws per component. For
multivariate: nested list with dimensions as names, each containing
variance draws for that dimension.

# Extract Covariance Draws from brms Model

Extracts posterior draws of covariances between dimensions for
multivariate models. Used for computing composite variance components.

## Usage

``` r
extract_covariance_draws(model, dimensions)
```

## Arguments

- model:

  A fitted brms model object

- dimensions:

  Character vector of dimension/response names

## Value

Named list with:

- residual:

  List of residual covariance draws (named by dimension pairs)

- random_effect:

  Named list of covariances per facet

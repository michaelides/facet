# Compute Observed Score Variances and Covariances for VAR Calculation

Computes observed score variances and covariances from the original
data, which are required for the Haberman (2008) formula in VAR
calculation. These are data quantities (not model estimates) and are
constant across posterior draws.

## Usage

``` r
compute_observed_covariances(data, dimensions, weights, dimension_var = NULL)
```

## Arguments

- data:

  A data frame containing the subscale scores.

- dimensions:

  Character vector of dimension/subscale names.

- weights:

  Named numeric vector of weights for each dimension.

## Value

A list with:

- obs_var:

  Named vector of observed score variances per dimension

- obs_cov:

  Matrix of observed score covariances between dimensions

- obs_var_C:

  Observed variance of the weighted composite

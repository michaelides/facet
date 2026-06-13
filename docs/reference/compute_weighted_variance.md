# Compute Weighted Variance from Covariance Matrix

Computes the variance of a weighted composite using the quadratic form
w' \* Sigma \* w, where w is the weight vector and Sigma is the
variance-covariance matrix.

## Usage

``` r
compute_weighted_variance(cov_matrix, weights)
```

## Arguments

- cov_matrix:

  A symmetric variance-covariance matrix

- weights:

  Numeric vector of weights (one per dimension)

## Value

The weighted variance as a numeric scalar

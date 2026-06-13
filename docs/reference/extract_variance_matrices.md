# Extract Variance-Covariance Matrices for Optimization

Extracts the universe score and error variance-covariance matrices from
a dstudy object for use in weight optimization.

## Usage

``` r
extract_variance_matrices(dstudy_obj, use_posterior_mean = TRUE)
```

## Arguments

- dstudy_obj:

  A dstudy object

- use_posterior_mean:

  If TRUE, use posterior means; if FALSE, use point estimates

## Value

List with Sigma_tau and Sigma_delta matrices

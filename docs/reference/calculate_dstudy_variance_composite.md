# Calculate D-Study Variance Components with Composite Row

Adds composite variance component rows for multivariate models using
posterior draws for full uncertainty propagation.

## Usage

``` r
calculate_dstudy_variance_composite(
  vc_draws,
  cov_draws,
  weights,
  n,
  object,
  n_provided
)
```

## Arguments

- vc_draws:

  Variance draws (from extract_variance_draws)

- cov_draws:

  Covariance draws (from extract_covariance_draws)

- weights:

  Named vector of weights

- n:

  Sample sizes (named list)

- object:

  Object specification

- n_provided:

  Whether n was provided

## Value

List with:

- vc_table:

  Tibble with variance components including composite rows

- posterior:

  List of composite draws

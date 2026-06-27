# Generate Variance and Covariance Draws from mom Point Estimates

Uses variance estimates and SEs to generate normal pseudo-posterior
draws, and converts correlation draws to covariance draws.

## Usage

``` r
generate_mom_variance_and_covariance_draws(gstudy_obj, n_draws = 1000)
```

## Arguments

- gstudy_obj:

  A gstudy object with mom estimator

- n_draws:

  Number of pseudo-draws to generate (default 1000)

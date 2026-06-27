# Extract Variance-Covariance Matrices from gstudy Objects

Extracts the variance-covariance matrices of random effects from a
gstudy object. This is a wrapper that calls the appropriate VarCorr
method based on the estimator used to fit the model.

## Usage

``` r
# S3 method for class 'gstudy'
VarCorr(x, ...)
```

## Arguments

- x:

  A gstudy object.

- ...:

  Additional arguments passed to the underlying VarCorr method (e.g.,
  [`lme4::VarCorr`](https://rdrr.io/pkg/nlme/man/VarCorr.html) or
  [`brms::VarCorr`](https://rdrr.io/pkg/nlme/man/VarCorr.html)).

## Value

For lme4 estimator, returns a list of variance-covariance matrices for
each random effect term. For brms estimator, returns a list with
posterior summaries.

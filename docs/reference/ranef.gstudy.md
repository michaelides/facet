# Extract Random Effects from gstudy Objects

Extracts the random effects (BLUPs/conditional modes) from a gstudy
object. This is a wrapper that calls the appropriate ranef method based
on the estimator used to fit the model.

## Usage

``` r
# S3 method for class 'gstudy'
ranef(object, ...)
```

## Arguments

- object:

  A gstudy object.

- ...:

  Additional arguments passed to the underlying ranef method (e.g.,
  [`lme4::ranef`](https://rdrr.io/pkg/lme4/man/ranef.html) or
  [`brms::ranef`](https://rdrr.io/pkg/nlme/man/random.effects.html)).

## Value

The random effects from the fitted model. For lme4 estimator, returns a
list of matrices. For brms estimator, returns a list of arrays with
posterior summaries.

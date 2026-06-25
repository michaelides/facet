# Define Priors Using Tidy Evaluation

This function is a wrapper around
[`brms::prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)
to allow specifying priors using tidy evaluation syntax. See
[`brms::prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)
for full documentation.

## Usage

``` r
prior(prior, ...)
```

## Arguments

- prior:

  A call or expression defining a distribution.

- ...:

  Additional arguments passed to
  [`brms::prior()`](https://paulbuerkner.com/brms/reference/set_prior.html).

## Value

A brmsprior object.

## References

Paul Bürkner (2017). brms: An R Package for Bayesian Multilevel Models.
*Journal of Statistical Software*, 80(1), 1-28.
[doi:10.18637/jss.v080.i01](https://doi.org/10.18637/jss.v080.i01)

## Examples

``` r
if (FALSE) { # \dontrun{
prior <- prior(normal(0, 1), class = sd)
g <- gstudy(score ~ (1 | person), data = mydata, prior = prior, backend = "brms",
  iter = 2000, cores = 4, refresh = 1000)
} # }
```

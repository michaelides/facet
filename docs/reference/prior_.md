# Define Priors Using Formula Syntax

This function is a wrapper around
[`brms::prior_()`](https://paulbuerkner.com/brms/reference/set_prior.html)
to allow specifying priors using formula syntax (with `~`). See
[`brms::prior_()`](https://paulbuerkner.com/brms/reference/set_prior.html)
for full documentation.

## Usage

``` r
prior_(prior, ...)
```

## Arguments

- prior:

  A one-sided formula defining the prior distribution.

- ...:

  Additional arguments passed to
  [`brms::prior_()`](https://paulbuerkner.com/brms/reference/set_prior.html).

## Value

A brmsprior object.

## References

Paul Bürkner (2017). brms: An R Package for Bayesian Multilevel Models.
*Journal of Statistical Software*, 80(1), 1-28.
[doi:10.18637/jss.v080.i01](https://doi.org/10.18637/jss.v080.i01)

## Examples

``` r
if (FALSE) { # \dontrun{
prior <- prior_(~ normal(0, 1), class = "sd")
g <- gstudy(score ~ (1 | person), data = mydata, prior = prior, backend = "brms",
  iter = 2000, cores = 4, refresh = 1000)
} # }
```

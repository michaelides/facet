# Get Default Priors for Bayesian Models

This function is a wrapper around
[`brms::default_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
to allow inspecting default priors that would be applied to a model. See
[`brms::default_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
for full documentation.

## Usage

``` r
default_prior(object, ...)
```

## Arguments

- object:

  An object whose class determines the method. Typically a formula,
  brmsformula, or mvbrmsformula.

- ...:

  Additional arguments passed to
  [`brms::default_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
  including `data` (a data frame containing the model data).

## Value

A brmsprior object; a data frame containing information about all
parameters for which priors can be specified.

## References

Paul Bürkner (2017). brms: An R Package for Bayesian Multilevel Models.
*Journal of Statistical Software*, 80(1), 1-28.
[doi:10.18637/jss.v080.i01](https://doi.org/10.18637/jss.v080.i01)

## Examples

``` r
if (FALSE) { # \dontrun{
# Get default priors for a gstudy model
prior_info <- default_prior(score ~ (1 | person) + (1 | item), data = mydata)
print(prior_info)
} # }
```

# Define Priors for brms Models

This function is a wrapper around
[`brms::set_prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)
to allow specifying priors for Bayesian models fit via the brms backend
in gstudy. See
[`brms::set_prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)
for full documentation on parameter specifications.

## Usage

``` r
set_prior(
  prior,
  class = "b",
  coef = "",
  group = "",
  resp = "",
  dpar = "",
  nlpar = "",
  lb = NA,
  ub = NA,
  check = TRUE,
  ...
)
```

## Arguments

- prior:

  A character string defining a distribution in Stan language.

- class:

  The parameter class. Defaults to `"b"` (i.e., population-level
  effects). See Details in
  [`brms::set_prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)
  for other valid parameter classes.

- coef:

  For class `"b"`, the coefficient name.

- group:

  For class `"sd"` or `"cor"`, the grouping factor.

- resp:

  For multivariate models, the response variable.

- dpar:

  For distributional parameters, the parameter name.

- nlpar:

  For non-linear parameters, the parameter name.

- lb:

  Lower bound for parameters with lower bounds.

- ub:

  Upper bound for parameters with upper bounds.

- check:

  Logical; check if prior is valid for the model.

- ...:

  Additional arguments passed to
  [`brms::set_prior()`](https://paulbuerkner.com/brms/reference/set_prior.html).

## Value

A brmsprior object; a data frame containing the prior specifications.

## References

Paul Bürkner (2017). brms: An R Package for Bayesian Multilevel Models.
*Journal of Statistical Software*, 80(1), 1-28.
[doi:10.18637/jss.v080.i01](https://doi.org/10.18637/jss.v080.i01)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define a prior for standard deviations of random effects
prior <- set_prior("normal(0, 1)", class = "sd", group = "person")

# Use with gstudy
g <- gstudy(score ~ (1 | person) + (1 | item),
  data = mydata,
  prior = prior,
  backend = "brms",
  iter = 2000, cores = 4, refresh = 1000
)
} # }
```

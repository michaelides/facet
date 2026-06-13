# Create an Empty Prior Object

This function is a wrapper around
[`brms::empty_prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)
to create an empty prior object to which priors can be added. See
[`brms::empty_prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)
for full documentation.

## Usage

``` r
empty_prior()
```

## Value

A brmsprior object with no rows.

## References

Paul Bürkner (2017). brms: An R Package for Bayesian Multilevel Models.
*Journal of Statistical Software*, 80(1), 1-28.
[doi:10.18637/jss.v080.i01](https://doi.org/10.18637/jss.v080.i01)

## Examples

``` r
if (FALSE) { # \dontrun{
# Start with empty prior and add custom priors
my_prior <- empty_prior()
my_prior <- rbind(my_prior, set_prior("normal(0, 1)", class = "sd", group = "person"))
g <- gstudy(score ~ (1 | person), data = mydata, prior = my_prior, backend = "brms")
} # }
```

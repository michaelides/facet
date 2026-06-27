# ANOVA-based Estimation Fit Object

An object representing a G-study model fitted using the method of
moments.

## Format

A list with class "aovfit" containing:

- formula:

  The model formula

- data:

  The data frame

- aov_model:

  The aov model object

- anova_results:

  List of ANOVA results per stratum

- variance_components:

  Tibble of estimated variance components

- response:

  Name of the response variable

- random_facets:

  Character vector of facet names

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, each = 5)),
  item = factor(rep(1:5, times = 20))
)
fit <- fit_aov(score ~ (1 | person) + (1 | item), data)
summary(fit)
} # }
```

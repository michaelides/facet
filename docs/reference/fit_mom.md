# Fit a Model Using Method of Moments

Fits a mixed effects model using the method of moments (ANOVA-based)
variance component estimation. Works best with balanced designs.

## Usage

``` r
fit_mom(formula, data, unbalanced = FALSE, nested = NULL, ...)
```

## Arguments

- formula:

  A formula for the model (lme4 style), e.g.,
  `y ~ (1 | facet1) + (1 | facet2)`.

- data:

  A data frame containing the variables in the formula.

- unbalanced:

  Logical. If TRUE, enables unbalanced multivariate estimation using
  Henderson's Method III. Each dimension is analyzed with its available
  data. Default is FALSE.

- ...:

  Additional arguments (currently unused).

## Value

An object of class "momfit" containing:

- formula:

  The formula used

- data:

  The data

- aov_model:

  The aov model object

- anova_results:

  The ANOVA results for each random effect

- variance_components:

  The estimated variance components as a tibble

- response:

  Name of the response variable

- random_facets:

  Character vector of random effect facet names

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, each = 5)),
  item = factor(rep(1:5, times = 20))
)
fit_mom(score ~ (1 | person) + (1 | item), data)
} # }
```

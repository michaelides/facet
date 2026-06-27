# Select the Appropriate Estimator

Determines which estimator to use based on formula characteristics and
user preference. If estimator is "auto", selects brms for multivariate
formulas and lme4 for univariate.

## Usage

``` r
select_estimator(formula, estimator = c("auto", "lme4", "brms", "aov"))
```

## Arguments

- formula:

  The model formula.

- estimator:

  Character string specifying the desired estimator: "auto", "lme4", or
  "brms".

## Value

Character string indicating the estimator to use.

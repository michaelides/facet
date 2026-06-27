# Validate Formula for G-Study

Checks that a formula is valid for a G-study analysis. The formula
should have a response variable and at least one random effect term.

## Usage

``` r
validate_formula(formula, estimator = "auto")
```

## Arguments

- formula:

  A formula object.

- estimator:

  Character string indicating the estimator to use ("auto", "lme4",
  "brms", or "aov").

## Value

TRUE if valid, otherwise raises an error.

# Convert Formula to Estimator-Specific Format

Converts a G-study formula to the format required by a specific
estimator. For lme4, this removes any brms-specific syntax.

## Usage

``` r
convert_formula(formula, estimator)
```

## Arguments

- formula:

  A formula object.

- estimator:

  Character string indicating the target estimator.

## Value

A formula object in the estimator-specific format.

# Detect lme4 Estimation Issues

Checks for convergence issues and singularity in lme4 models.

## Usage

``` r
detect_lme4_issues(model)
```

## Arguments

- model:

  An lmerMod object from lme4::lmer().

## Value

A list with components:

- convergence:

  Logical indicating if convergence issues were detected

- singularity:

  Logical indicating if the model is singular
